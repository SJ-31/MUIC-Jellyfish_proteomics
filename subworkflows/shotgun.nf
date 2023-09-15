// CONVENTION: the final protein/PSMs output file will be prefixed with the engine that produced them
include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { COMET } from '../modules/comet'
include { IDENTIPY; FORMAT_IDPY } from '../modules/identipy'
include { METAMORPHEUS as METAMORPHEUS_DEFAULT } from '../modules/metamorpheus'
include { METAMORPHEUS as METAMORPHEUS_GLYCO } from '../modules/metamorpheus'
include { TIDE } from '../modules/tide'
include { TIDE_COMBINED_PEP } from '../modules/tide'
include { FORMAT_MQ } from '../modules/format_mq'
include { PERCOLATOR } from '../modules/percolator'
include { DEISOTOPE } from '../modules/deisotope'
include { MS_MAPPING } from '../modules/ms_mapping'
include { CALIBRATE } from '../modules/calibrate'
include { FALCON } from '../modules/falcon'
include { bk_decoys } from './bk_decoys.nf'
include { combine_searches as combine_searches_FIRST } from './combine_searches'
include { combine_searches as combine_searches_SECOND } from './combine_searches'
include { quantify as quantify_FIRST } from './quantify'
include { quantify as quantify_SECOND } from './quantify'
include { make_db } from './make_db'
include { open_search as open_search_FIRST } from './open_search'
include { open_search as open_search_SECOND } from './open_search'

workflow 'pre' {
    take:
    mzML
    raw
    normal_database

    main:
    DEISOTOPE(mzML.collect(),"$params.results/Preprocessed")
    FALCON(mzML.collect(), "$params.results/Preprocessed/Falcon")
    CALIBRATE(raw.collect(), normal_database,
              "$params.results/Preprocessed/Falcon")
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
    db_normal

    main:
    Channel.fromPath(params.db_spec).splitText() { it.replaceAll("\n", "") }
        .branch {
            normal: it ==~ /.*all_normal.fasta/
            plusdecoys: it ==~ /.*decoysWnormal.fasta/
            seq_header_mapping : it ==~ /.*mappings.tsv/
            downloaded : it ==~ /.*downloaded.fasta/
        }.set { db }

    MS_MAPPING(mzML.collect(), "$params.results")

    // All searches
    MAXQUANT(raw.collect(), "default_maxquant.xml",
             "$params.results/1-First_pass/Engines/MaxQuant",
             "$params.results/1-First_pass/Logs", db.normal.first())
    COMET(mzML.collect(), "$params.results/1-First_pass/Engines/Comet",
          "$params.results/1-First_pass/Logs", db.plusdecoys)
    MSFRAGGER(mzML.collect(), "$params.config/MSFragger_params.params", "",
              "$params.results/1-First_pass/Engines/MsFragger",
              "$params.results/1-First_pass/Logs", db.plusdecoys)
    IDENTIPY(mzML.collect(), "$params.results/1-First_pass/Engines/Identipy",
             "$params.results/1-First_pass/Logs", db.plusdecoys.first())
    FORMAT_IDPY(IDENTIPY.out.pepxml.collect(),
                "$params.results/1-First_pass/Engines/Identipy")
    METAMORPHEUS_DEFAULT(mgf.collect(),
                         "$params.results/1-First_pass/Engines/Metamorpheus",
                         "$params.results/1-First_pass/Logs", "",
                         "$params.config/metamorpheus_params.toml", db.normal)
    // MSGF(mzML, "$params.results/1-First_pass/msgf", db.normal) No longer used,
    //  no way to integrate with percolator for now
    TIDE(mgf.collect(), "$params.results/1-First_pass/Engines/Tide",
         "$params.results/1-First_pass/Logs",
         "$params.results/1-First_pass/Percolator",
    db.normal)
    TIDE_COMBINED_PEP(TIDE.out.percolator, "$params.results/1-First_pass/Percolator")
    FORMAT_MQ(MAXQUANT.out.msmsScans.collect(), "$params.results/1-First_pass/Engines/MaxQuant")

    // Post-processing with Percolator
    empty = Channel.empty()
    empty.mix(
        METAMORPHEUS_DEFAULT.out.percolator,
        FORMAT_MQ.out,
        COMET.out.percolator,
        MSFRAGGER.out.percolator,
        FORMAT_IDPY.out.percolator
    ).set { to_percolator }
    PERCOLATOR(to_percolator, "$params.results/1-First_pass/Percolator",
               "$params.results/1-First_pass/Logs",
               db.plusdecoys.first())

    // First pass quantification
    quantify_FIRST(MS_MAPPING.out,
                   mzML,
                   PERCOLATOR.out.psms.mix(TIDE.out.perc_psms),
                   METAMORPHEUS_DEFAULT.out.psms,
                   TIDE.out.target,
                   PERCOLATOR.out.psm2combinedPEP
                       .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP),
                   "$params.results/1-First_pass/Quantify")

    // First open search
    open_search_FIRST(quantify_FIRST.out.unmatched_msms,
                      db.plusdecoys,
                      db.normal,
                      "$params.results/1-First_pass/Open_search")

    // First combining
    combine_searches_FIRST(
        PERCOLATOR.out.prot2intersect
            .mix(TIDE.out.perc_protein).collect(),
        PERCOLATOR.out.psm2combinedPEP
            .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP).collect(),
        quantify_FIRST.out.directlfq,
        PERCOLATOR.out.psms.mix(TIDE.out.perc_psms),
        "$params.results/1-First_pass",
        db.seq_header_mapping)

    // Second pass with Bern and Kil decoy database
    compatible = /.*comet.*|.*identipy.*|.*msfragger.*/
    bk_decoys(PERCOLATOR.out.prot.filter({ it[0] =~ compatible }),
              db.seq_header_mapping, mzML)
    from_first = /.*metamorpheus.*|.*maxquant.*/

    // Second pass quantification
    quantify_SECOND(MS_MAPPING.out,
                    mzML,
                    bk_decoys.out.all_psms.mix(PERCOLATOR.out.psms
                                               .filter( ~from_first ),
                                               TIDE.out.perc_psms),
                   METAMORPHEUS_DEFAULT.out.psms,
                   TIDE.out.target,
                   bk_decoys.out.psm2combinedPEP.mix(
                       PERCOLATOR.out.psm2combinedPEP.filter( ~from_first ),
                       TIDE_COMBINED_PEP.out.psm2combinedPEP),
                   "$params.results/2-Second_pass/Quantify")

    // Second open search
    open_search_SECOND(quantify_SECOND.out.unmatched_msms,
                      db.plusdecoys,
                      db.normal,
                      "$params.results/1-Second_pass/Open_search")

    // Second combining
    combine_searches_SECOND(
        bk_decoys.out.prot2intersect.mix(
            PERCOLATOR.out.prot2intersect.filter( ~from_first ),
            TIDE.out.perc_protein).collect(),
        bk_decoys.out.psm2combinedPEP.mix(
            PERCOLATOR.out.psm2combinedPEP.filter( ~from_first ),
            TIDE_COMBINED_PEP.out.psm2combinedPEP).collect(),
       quantify_SECOND.out.directlfq,
       bk_decoys.out.all_psms.mix(PERCOLATOR.out.psms
                                  .filter( ~from_first ),
                                  TIDE.out.perc_psms),
       "$params.results/2-Second_pass",
       db.seq_header_mapping)
}
