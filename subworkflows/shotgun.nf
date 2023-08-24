// CONVENTION: the final protein/PSMs output file will be prefixed with the engine that produced them
include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { COMET } from '../modules/comet'
include { IDENTIPY } from '../modules/identipy'
include { METAMORPHEUS } from '../modules/metamorpheus'
include { TIDE } from '../modules/tide'
include { TIDE_COMBINED_PEP } from '../modules/tide'
include { PERCOLATOR } from '../modules/percolator'
include { MS2RESCORE } from '../modules/ms2rescore'
include { DEISOTOPE } from '../modules/deisotope'
include { MS_MAPPING } from '../modules/ms_mapping'
include { QUANDENSER } from '../modules/quandenser'
include { FALCON } from '../modules/falcon'
include { bk_decoys } from './bk_decoys.nf'
include { combine_searches as combine_searches_FIRST } from './combine_searches'
include { combine_searches as combine_searches_SECOND } from './combine_searches'
include { quantify as quantify_FIRST } from './quantify'
include { quantify as quantify_SECOND } from './quantify'
include { make_db } from './make_db'


Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.indexed_mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mgf: it =~ /.mgf/
        mzXML: it =~ /.mzXML/
        raw: it =~ /.raw/
    }.set { manifest }


workflow 'preprocess' {
    manifest.mzML.collectFile("$params.results/Preprocessed/manifest.txt",
                              newLine: true)
        .set { mzML_spec }
    DEISOTOPE(mzML_spec,
              "$params.results/Preprocessed/Proteowizard")
    QUANDENSER(mzML_spec,
              "$params.results/Preprocessed/Quandenser")
    FALCON(manifest.mzML, "$params.results/Preprocessed/Falcon")
}


workflow 'search' {

    take:
    db_list

    main:
    if ( params.with_db ) {
        databases = make_db().out
    } else {
        databases = params.db_spec
    }
    Channel.fromPath(databases).splitText() { it.replaceAll("\n", "") }
        .branch {
            normal: it ==~ /.*all_normal.fasta/
            plusdecoys: it ==~ /.*decoysWnormal.fasta/
            seq_mapping : it ==~ /.*decoysWnormal_mapping.tsv/
            header_mapping : it ==~ /.*header_mappings.tsv/
        }.set { db }
    empty = Channel.empty()

    MS_MAPPING(manifest.mzML.collect(), "$params.results")

    // All searches
    MAXQUANT(manifest.raw, "$params.results/1-First_pass/MaxQuant", db.normal.first())
    COMET(manifest.mzML.collect(), "$params.results/1-First_pass/Comet",
    db.plusdecoys)
    MSFRAGGER(manifest.mzML.collect(), "$params.config/MSFragger_params.params",
    "$params.results/1-First_pass/MsFragger", db.plusdecoys)
    IDENTIPY(manifest.mzML.collect(), "$params.results/1-First_pass/Identipy", db.plusdecoys)
    METAMORPHEUS(manifest.mgf.collect(), "$params.results/1-First_pass/Metamorpheus", db.normal)
    // MSGF(manifest.mzML, "$params.results/1-First_pass/msgf", db.normal) No longer used,
    //  no way to integrate with percolator for now
    TIDE(manifest.mgf.collect(), "$params.results/1-First_pass/Tide", "$params.results/1-First_pass/Percolator",
    db.normal)
    TIDE_COMBINED_PEP(TIDE.out.percolator, "$params.results/1-First_pass/Percolator")
    MAXQUANT.out.ms2rescore.collect()
        .map { it -> ["maxquant", it] }
        .set { maxqms2rescore }
    MS2RESCORE(maxqms2rescore, "$params.results/1-First_pass/MaxQuant",
    manifest.mgf.collect())

    // Post-processing with Percolator
    empty.mix( // Perhaps you need to start with ms2rescore
        METAMORPHEUS.out.percolator,
        MS2RESCORE.out.pin,
        COMET.out.percolator,
        MSFRAGGER.out.percolator,
        IDENTIPY.out.percolator
    ).set { to_percolator }
    PERCOLATOR(to_percolator, "$params.results/1-First_pass/Percolator", db.plusdecoys.first())

    quantify_FIRST(MS_MAPPING.out, manifest.mzML,
                   PERCOLATOR.out.psms.mix(TIDE.out.perc_psms),
                   METAMORPHEUS.out.psms,
                   MS2RESCORE.out.pin.flatten().filter(~/.*pin.*/),
                   TIDE.out.target,
                   "$params.results/1-First_pass/1-Quantify")
    // First combining
    combine_searches_FIRST(
        PERCOLATOR.out.prot2intersect
            .mix(TIDE.out.perc_protein).collect(),
        PERCOLATOR.out.psm2combinedPEP
            .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP).collect(),
        PERCOLATOR.out.prot2combinedPEP
            .mix(TIDE_COMBINED_PEP.out.prot2combinedPEP).collect(),
        "$params.results/1-First_pass",
        db.header_mapping,
        db.seq_mapping)

    // Second pass with Bern and Kil decoy database
    compatible = /.*comet.*|.*identipy.*|.*msfragger.*/
    bk_decoys(PERCOLATOR.out.prot.filter({ it[0] =~ compatible }),
              db.seq_mapping, db.header_mapping, manifest.mzML)
    from_first = /.*metamorpheus.*|.*maxquant.*/

    quantify_SECOND(MS_MAPPING.out, manifest.mzML,
                    bk_decoys.out.all_psms.mix(PERCOLATOR.out.psms
                                               .filter( ~from_first ),
                                               TIDE.out.perc_psms),
                   METAMORPHEUS.out.psms,
                    MS2RESCORE.out.pin.flatten().filter(~/.*pin.*/),
                   TIDE.out.target,
                   "$params.results/2-Second_pass/1-Quantify")

    combine_searches_SECOND(
        bk_decoys.out.prot2intersect.mix(
            PERCOLATOR.out.prot2intersect.filter( ~from_first ),
            TIDE.out.perc_protein).collect(),
        bk_decoys.out.psm2combinedPEP.mix(
            PERCOLATOR.out.psm2combinedPEP.filter( ~from_first ),
            TIDE_COMBINED_PEP.out.psm2combinedPEP).collect(),
        bk_decoys.out.prot2combinedPEP.mix(
            PERCOLATOR.out.prot2combinedPEP.filter( ~from_first ),
            TIDE_COMBINED_PEP.out.prot2combinedPEP).collect(),
        "$params.results/2-Second_pass",
        db.header_mapping,
        db.seq_mapping)
}
