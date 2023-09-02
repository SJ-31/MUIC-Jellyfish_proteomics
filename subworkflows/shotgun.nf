// CONVENTION: the final protein/PSMs output file will be prefixed with the engine that produced them
include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { COMET } from '../modules/comet'
include { IDENTIPY } from '../modules/identipy'
include { METAMORPHEUS } from '../modules/metamorpheus'
include { TIDE } from '../modules/tide'
include { TIDE_COMBINED_PEP } from '../modules/tide'
include { FORMAT_MQ } from '../modules/format_mq'
include { PERCOLATOR } from '../modules/percolator'
include { DEISOTOPE } from '../modules/deisotope'
include { MS_MAPPING } from '../modules/ms_mapping'
include { FALCON } from '../modules/falcon'
include { bk_decoys } from './bk_decoys.nf'
include { combine_searches as combine_searches_FIRST } from './combine_searches'
include { combine_searches as combine_searches_SECOND } from './combine_searches'
include { quantify as quantify_FIRST } from './quantify'
include { quantify as quantify_SECOND } from './quantify'
include { make_db } from './make_db'

workflow 'pre' {
    take:
    mzML

    main:
    DEISOTOPE(mzML.collect(),"$params.results/Preprocessed")
    FALCON(mzML.collect(), "$params.results/Preprocessed/Falcon")
}

println """
Searching spectra...
    Prefix: $params.pref
    Manifest file: $params.manifest_file
    Database file: $params.db_spec
    Results path: $params.results
"""

workflow 'search' {

    take:
    mzML
    mgf
    raw
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
            downloaded : it ==~ /.*downloaded.fasta/
        }.set { db }
    empty = Channel.empty()

    MS_MAPPING(mzML.collect(), "$params.results")

    // All searches
    MAXQUANT(raw, "$params.results/1-First_pass/MaxQuant", db.normal.first())
    COMET(mzML.collect(), "$params.results/1-First_pass/Engines/Comet",
    db.plusdecoys)
    MSFRAGGER(mzML.collect(), "$params.config/MSFragger_params.params",
    "$params.results/1-First_pass/Engines/MsFragger", db.plusdecoys)
    IDENTIPY(mzML.collect(), "$params.results/1-First_pass/Identipy", db.plusdecoys)
    METAMORPHEUS(mgf.collect(), "$params.results/1-First_pass/Engines/Metamorpheus", db.normal)
    // MSGF(mzML, "$params.results/1-First_pass/msgf", db.normal) No longer used,
    //  no way to integrate with percolator for now
    TIDE(mgf.collect(), "$params.results/1-First_pass/Engines/Tide", "$params.results/1-First_pass/Percolator",
    db.normal)
    TIDE_COMBINED_PEP(TIDE.out.percolator, "$params.results/1-First_pass/Percolator")
    FORMAT_MQ(MAXQUANT.out.msmsScans.collect(), "$params.results/1-First_pass/Engines/MaxQuant")

    // Post-processing with Percolator
    empty.mix(
        METAMORPHEUS.out.percolator,
        FORMAT_MQ.out,
        COMET.out.percolator,
        MSFRAGGER.out.percolator,
        IDENTIPY.out.percolator
    ).set { to_percolator }
    PERCOLATOR(to_percolator, "$params.results/1-First_pass/Percolator", db.plusdecoys.first())

    quantify_FIRST(MS_MAPPING.out, mzML,
                   PERCOLATOR.out.psms.mix(TIDE.out.perc_psms),
                   METAMORPHEUS.out.psms,
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
        quantify_FIRST.out.directlfq,
        "$params.results/1-First_pass",
        db.header_mapping,
        db.seq_mapping)

    // Second pass with Bern and Kil decoy database
    compatible = /.*comet.*|.*identipy.*|.*msfragger.*/
    bk_decoys(PERCOLATOR.out.prot.filter({ it[0] =~ compatible }),
              db.seq_mapping, db.header_mapping, mzML)
    from_first = /.*metamorpheus.*|.*maxquant.*/

    quantify_SECOND(MS_MAPPING.out, mzML,
                    bk_decoys.out.all_psms.mix(PERCOLATOR.out.psms
                                               .filter( ~from_first ),
                                               TIDE.out.perc_psms),
                   METAMORPHEUS.out.psms,
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
        quantify_SECOND.out.directlfq,
        "$params.results/2-Second_pass",
        db.header_mapping,
        db.seq_mapping)
}
