// CONVENTION: the final protein/PSMs output file will be prefixed with the engine that produced them
include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { COMET } from '../modules/comet'
include { IDENTIPY } from '../modules/identipy'
include { METAMORPHEUS } from '../modules/metamorpheus'
include { TIDE } from '../modules/tide'
include { TIDE_COMBINED_PEP } from '../modules/tide'
include { CASANOVO } from '../modules/casanovo'
include { PERCOLATOR } from '../modules/percolator'
include { MS2RESCORE } from '../modules/ms2rescore'
include { COMBINED_DATABASE } from '../modules/combined_database'
include { EXTRACT_CASANOVO } from '../modules/casanovo'
include { DEISOTOPE } from '../modules/deisotope'
include { PEPNET } from '../modules/pepnet'
include { EXTRACT_PEPNET } from '../modules/pepnet'
include { bk_decoys } from './bk_decoys.nf'
include { combine_searches as combine_searches_FIRST } from './combine_searches.nf'
include { combine_searches as combine_searches_SECOND } from './combine_searches.nf'

workflow 'preprocess' {
    DEISOTOPE(params.manifest_file)
}

Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mzXML: it =~ /.mzXML/
        mgf: it =~ /.mgf/
        raw: it =~ /.raw/
    }.set { manifest }

Channel.fromPath(params.databases)
    .splitText()
    .set { database_listing }

workflow 'make_db' {
    if ( params.denovo ) {
    // SMSNET(manifest.mgf.collect(), "$params.results/SMSNET") // TODO: Fixthis
    CASANOVO(manifest.mzML,"$params.results/Casanovo")
    EXTRACT_CASANOVO(CASANOVO.out.collect(), "$params.results/Casanovo")
    PEPNET(manifest.mgf, "$params.results/PepNet")
    EXTRACT_PEPNET(PEPNET.out.collect(), "$params.results/PepNet")

    EXTRACT_CASANOVO.out.mix(EXTRACT_PEPNET.out)
     .set { denovo }
    } else {
        denovo = Channel.empty()
    }
    COMBINED_DATABASE(database_listing.collect(), denovo.collect(),
                      "$params.results/databases")
}

workflow 'search' {
    Channel.fromPath(params.db_spec)
        .splitText()
        .flatten().branch {
            normal: it =~ /all_normal/
            plusdecoys: it =~ /decoysWnormal.fasta/
            seq_mapping : it =~ /decoysWnormal_mapping/
            header_mapping : it =~ /header_mappings/
        }.set { db }
    empty = Channel.empty()
    // MaxQuant seems to only work with .raw files
    MAXQUANT(manifest.raw, "$params.results/MaxQuant", db.normal)
    COMET(manifest.mzXML.collect(), "$params.results/First_pass/Comet",
    db.plusdecoys)
    MSFRAGGER(manifest.mzML.collect(), "$params.config/MSFragger_params.params",
    "$params.results/First_pass/MsFragger", db.plusdecoys)
    IDENTIPY(manifest.mzML.collect(), "$params.results/First_pass/Identipy", db.plusdecoys)
    METAMORPHEUS(manifest.mzML.collect(), "$params.results/First_pass/Metamorpheus", db.normal)
    // MSGF(manifest.mzML, "$params.results/First_pass/msgf", db.normal) No longer used,
    //  no way to integrate with percolator for now
    TIDE(manifest.mzXML.collect(), "$params.results/First_pass/Tide", "$params.results/First_pass/Percolator",
    db.normal)
    TIDE_COMBINED_PEP(TIDE.out.percolator, "$params.results/First_pass/Percolator")
    TIDE.out.percolator.flatten().filter( ~/.*_target_proteins\.tsv/ )
        .set { tide_percolator }
    MAXQUANT.out.ms2rescore.flatten().filter( ~/.*txt/ ).collect()
        .map { it -> ["maxquant", it] }
        .set { maxqms2rescore }
    MS2RESCORE(maxqms2rescore, "$params.results/First_pass/MaxQuant",
    manifest.mgf.collect())
    empty.mix(
        METAMORPHEUS.out.percolator,
        MS2RESCORE.out,
        COMET.out.percolator,
        MSFRAGGER.out.percolator,
        IDENTIPY.out.percolator
    ).set { to_percolator }
    PERCOLATOR(to_percolator, "$params.results/First_pass/Percolator", db.plusdecoys)
    combine_searches_FIRST(PERCOLATOR.out.prot2intersect.mix(tide_percolator),
                           PERCOLATOR.out.psm2combinedPEP
                                .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP),
                           PERCOLATOR.out.prot2combinedPEP
                                .mix(TIDE_COMBINED_PEP.out.prot2combinedPEP),
                        "$params.results/First_pass")

    // Second pass with Bern and Kil decoy database
    compatible = /.*comet.*|.*identipy.*|.*msfragger.*/
    bk_decoys(PERCOLATOR.out.prot.filter( ~compatible), db.seq_mapping,
              db.header_mapping, manifest.mzXML, manifest.mzML)
    from_first = /.*metamorpheus.*|.*maxquant.*/
    combine_searches_SECOND(
        PERCOLATOR.out.prot2intersect.filter( ~from_first )
            .mix(tide_percolator, bk_decoys.out.prot2intersect),
        PERCOLATOR.out.psm2combinedPEP.filter( ~from_first )
            .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP,
                 bk_decoys.out.psm2combinedPEP),
        PERCOLATOR.out.prot2combinedPEP.filter( ~from_first )
            .mix(TIDE_COMBINED_PEP.out.prot2combinedPEP,
                bk_decoys.out.prot2combinedPEP),
        db.header_mapping,
        "$params.results/Second_pass")
}
