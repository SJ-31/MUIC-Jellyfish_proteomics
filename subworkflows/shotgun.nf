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
    .map { it -> [ it.Prefix, it.Raw, it.indexed_mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mgf: it =~ /.mgf/
        mzXML: it =~ /.mzXML/
        raw: it =~ /.raw/
    }.set { manifest }

workflow 'make_db' {

    main:
    Channel.fromPath(params.to_construct)
        .splitText()
        .set { database_listing }
    if ( params.denovo ) {
    // SMSNET(manifest.mgf.collect(), "$params.results/SMSNET") // TODO: Fixthis
    CASANOVO(manifest.mzML,"$params.results/Casanovo")
        EXTRACT_CASANOVO(CASANOVO.out.peps.collect(), "$params.results/Casanovo")
    PEPNET(manifest.mgf, "$params.results/PepNet")
    EXTRACT_PEPNET(PEPNET.out.peps.collect(), "$params.results/PepNet")

    EXTRACT_CASANOVO.out.mix(EXTRACT_PEPNET.out)
     .set { denovo }
    } else {
        denovo = Channel.empty()
    }
    COMBINED_DATABASE(database_listing.collect(), denovo.collect(),
                      "$params.results/databases")

    emit:
    COMBINED_DATABASE.out.listing
}

workflow 'search' {

    take:
    db_list

    main:
    if ( params.with_db ) {
        databases = db_list
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
    db.normal.view()
    empty = Channel.empty()

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
    TIDE.out.percolator.flatten().filter( ~/.*_percolator_proteins\.tsv/ )
        .set { tide_percolator }
    MAXQUANT.out.ms2rescore.flatten().filter( ~/.*txt/ ).collect()
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

    // First combining
    combine_searches_FIRST(
        PERCOLATOR.out.prot2intersect
            .mix(tide_percolator).collect(),
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
    combine_searches_SECOND(
        bk_decoys.out.prot2intersect.mix(
            PERCOLATOR.out.prot2intersect.filter( ~from_first ),
            tide_percolator).collect(),
        bk_decoys.out.psm2combinedPEP.mix(
            PERCOLATOR.out.psm2combinedPEP.filter( ~from_first ),
            TIDE_COMBINED_PEP.out.psm2combinedPEP).collect(),
        bk_decoys.out.prot2combinedPEP.mix(
            PERCOLATOR.out.prot2combinedPEP.filter( ~from_first ),
            TIDE_COMBINED_PEP.out.prot2combinedPEP).collect(),
        "$params.results/2-Second_pass",
        db.header_mapping,
        db.seq_mapping
        )
}
