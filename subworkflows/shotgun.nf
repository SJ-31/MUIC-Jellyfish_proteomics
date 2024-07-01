// CONVENTION: the final protein/PSMs output file will be prefixed with the engine that produced them
preprocessed_logs = "$params.results/Preprocessed/Logs"

include { MSGF; MSGF_MERGE } from  '../modules/msgf'
include { MSFRAGGER } from '../modules/msfragger'
include { COMET } from '../modules/comet'
include { IDENTIPY; FORMAT_IDPY } from '../modules/identipy'
include { METAMORPHEUS as METAMORPHEUS_DEFAULT } from '../modules/metamorpheus'
include { TIDE } from '../modules/tide'
include { TIDE_COMBINED_PEP } from '../modules/tide'
include { PERCOLATOR } from '../modules/percolator'
include { RAWPARSE } from '../modules/rawparse' addParams(logdir: preprocessed_logs)
include { DEISOTOPE } from '../modules/deisotope' addParams(logdir: preprocessed_logs)
include { MS_MAPPING } from '../modules/ms_mapping' addParams(logdir: preprocessed_logs)
include { CALIBRATE } from '../modules/calibrate' addParams(logdir: preprocessed_logs)
include { FALCON } from '../modules/falcon' addParams(logdir: preprocessed_logs)
include { bk_decoys } from './bk_decoys.nf'
include { combine_searches as combine_searches_FIRST } from './combine_searches'
include { combine_searches as combine_searches_SECOND } from './combine_searches'
include { quantify as quantify_FIRST } from './quantify'
include { quantify as quantify_SECOND } from './quantify'
include { open_search as open_search_FIRST } from './open_search'
include { open_search as open_search_SECOND } from './open_search'
include { FILTER_DB } from '../modules/filter_db'

workflow 'pre' {
    take:
    raw
    mzML
    normal_database

    main:
    outdir = "$params.results/Preprocessed"

    RAWPARSE(raw.collect(),"$outdir/Converted")
    DEISOTOPE(RAWPARSE.out.collect(),"$outdir/Proteowizard")
    FALCON(RAWPARSE.out.collect(), "$outdir/Falcon")
    CALIBRATE(mzML.collect(), normal_database, "$outdir/Metamorpheus")
}

workflow 'search' {

    take:
    mzML
    mgf
    db_normal
    db_plusdecoys
    db_seq_header_mapping


    main:
    MS_MAPPING(mzML.collect(), "$params.results")

    // All searches
    COMET(mzML.collect(), "$params.p1/Engines/Comet", "$params.p1/Logs",
          db_plusdecoys)
    MSFRAGGER(mzML.collect(), "$params.config_dir/MSFragger_params.params", "",
              "$params.p1/Engines/MsFragger", "$params.p1/Logs", db_plusdecoys)
    IDENTIPY(mzML.collect(), "$params.p1/Engines/Identipy", "$params.p1/Logs",
             db_plusdecoys.first())
    FORMAT_IDPY(IDENTIPY.out.pepxml.collect(),
                "$params.p1/Engines/Identipy")
    METAMORPHEUS_DEFAULT(mgf.collect(), "$params.p1/Engines/Metamorpheus",
                         "$params.p1/Logs", "",
                         "$params.config_dir/metamorpheus_params.toml", db_normal)
    MSGF(mzML, "$params.p1/Engines/MSGF", "$params.p1/Logs",
         db_plusdecoys.first())
    MSGF_MERGE(MSGF.out.pin.collect(), "$params.p1/Engines/MSGF")
    TIDE(mgf.collect(), "$params.p1/Engines/Tide", "$params.p1/Logs",
         "$params.p1/Percolator",
         db_normal)
    TIDE_COMBINED_PEP(TIDE.out.percolator, "$params.p1/Percolator")

    // Post-processing with Percolator
    empty = Channel.empty()
    empty.mix(
        METAMORPHEUS_DEFAULT.out.percolator,
        MSGF_MERGE.out,
        COMET.out.percolator,
        MSFRAGGER.out.percolator,
        FORMAT_IDPY.out.percolator
    ).set { to_percolator }
    PERCOLATOR(to_percolator, "$params.p1/Percolator", "$params.p1/Logs",
               db_plusdecoys.first())

    quantify_FIRST(MS_MAPPING.out,
                   mzML,
                   PERCOLATOR.out.psms.mix(TIDE.out.perc_psms),
                   PERCOLATOR.out.prot2intersect.mix(TIDE.out.perc_protein),
                   METAMORPHEUS_DEFAULT.out.psms,
                   TIDE.out.target,
                   PERCOLATOR.out.psm2combinedPEP
                       .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP),
                   "$params.p1/Quantify")

    open_search_FIRST(quantify_FIRST.out.unmatched_msms,
                      db_plusdecoys,
                      db_normal,
                      "$params.p1/Open_search",
                      db_seq_header_mapping)

    combine_searches_FIRST(
        PERCOLATOR.out.prot2intersect
            .mix(TIDE.out.perc_protein).collect(),
        PERCOLATOR.out.psm2combinedPEP
            .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP).collect(),
        quantify_FIRST.out.flashlfq,
        quantify_FIRST.out.maxlfq,
        quantify_FIRST.out.directlfq,
        quantify_FIRST.out.unmatched_pep_tsv,
        "$params.p1",
        db_seq_header_mapping,
        open_search_FIRST.out.open_results,
        open_search_FIRST.out.prot2intersect.collect(),
        quantify_FIRST.out.directlfq_input,
        params.blast_db)

    /* Second pass with Bern and Kil decoy database
    */
    compatible = /.*comet.*|.*identipy.*|.*msfragger.*|.*msgf.*/
    bk_decoys(PERCOLATOR.out.prot.filter({ it[0] =~ compatible }),
              db_seq_header_mapping, mzML)
    from_first = /.*metamorpheus.*/

    quantify_SECOND(MS_MAPPING.out,
                   mzML,
                   bk_decoys.out.all_psms.mix(PERCOLATOR.out.psms
                                              .filter( ~from_first ),
                                              TIDE.out.perc_psms),
                   bk_decoys.out.prot2intersect.mix(PERCOLATOR.out.prot2intersect
                                               .filter( ~from_first ),
                                               TIDE.out.perc_protein),
                   METAMORPHEUS_DEFAULT.out.psms,
                   TIDE.out.target,
                   bk_decoys.out.psm2combinedPEP.mix(
                       PERCOLATOR.out.psm2combinedPEP.filter( ~from_first ),
                       TIDE_COMBINED_PEP.out.psm2combinedPEP),
                   "$params.p2/Quantify")

    open_search_SECOND(quantify_SECOND.out.unmatched_msms,
                      db_plusdecoys,
                      db_normal,
                      "$params.p2/Open_search",
                      db_seq_header_mapping)

    combine_searches_SECOND(
        bk_decoys.out.prot2intersect.mix(
            PERCOLATOR.out.prot2intersect.filter( ~from_first ),
            TIDE.out.perc_protein).collect(),
        bk_decoys.out.psm2combinedPEP.mix(
            PERCOLATOR.out.psm2combinedPEP.filter( ~from_first ),
            TIDE_COMBINED_PEP.out.psm2combinedPEP).collect(),
       quantify_SECOND.out.flashlfq,
       quantify_SECOND.out.maxlfq,
       quantify_SECOND.out.directlfq,
       quantify_SECOND.out.unmatched_pep_tsv,
       "$params.p2",
       db_seq_header_mapping,
       open_search_SECOND.out.open_results,
       open_search_SECOND.out.prot2intersect.collect(),
       quantify_SECOND.out.directlfq_input,
       params.blast_db)
}
