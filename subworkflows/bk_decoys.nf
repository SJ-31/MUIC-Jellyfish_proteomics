logs = "$params.p2/Logs"

include { COMET } from  '../modules/comet'
include { MSGF; MSGF_MERGE } from  '../modules/msgf'
include { MSFRAGGER } from '../modules/msfragger'
include { PERCOLATOR as P_COMET } from '../modules/percolator'
include { PERCOLATOR as P_IPY} from '../modules/percolator'
include { PERCOLATOR as P_FRAGGER } from '../modules/percolator'
include { PERCOLATOR as P_MSGF } from '../modules/percolator'
include { IDENTIPY; FORMAT_IDPY } from '../modules/identipy'
include { MAKE_BK_DB } from '../modules/make_bk_db'

workflow bk_decoys {
    take:
    percolator_out
    seq_header_mapping
    mzML_ch

    main:
    MAKE_BK_DB(percolator_out, seq_header_mapping.first(),
               "$params.results/2-Second_pass/BK_databases")
    .flatten().filter( ~/.*\.fasta/ )
    .branch {
        comet: it.baseName =~ /comet*/
        identipy: it.baseName =~ /identipy*/
        msfragger: it.baseName =~ /msfragger*/
        msgf: it.baseName =~ /msgf*/
        }.set { bk_db }
    COMET(mzML_ch.collect(), "$params.p2/Engines/Comet", logs, bk_db.comet)
    IDENTIPY(mzML_ch.collect(), "$params.p2/Engines/Identipy", logs,
             bk_db.identipy.first())
    FORMAT_IDPY(IDENTIPY.out.pepxml,
                "$params.p2/Engines/Identipy")
    MSFRAGGER(mzML_ch.collect(), "$params.config_dir/MSFragger_params.params", "",
              "$params.p2/Engines/MsFragger", logs, bk_db.msfragger)
    MSGF(mzML_ch, "$params.p2/Engines/MSGF", logs, bk_db.msgf.first())
    MSGF_MERGE(MSGF.out.pin.collect(), "$params.p1/Engines/MSGF")
    P_MSGF(MSGF_MERGE.out,
            "$params.p2/Percolator", logs, bk_db.msgf)
    P_COMET(COMET.out.percolator,
            "$params.p2/Percolator", logs, bk_db.comet)
    P_IPY(FORMAT_IDPY.out.percolator,
          "$params.p2/Percolator", logs, bk_db.identipy)
    P_FRAGGER(MSFRAGGER.out.percolator,
              "$params.p2/Percolator", logs, bk_db.msfragger)

    emit:
    all_psms = P_COMET.out.psms
        .mix(P_IPY.out.psms,
             P_MSGF.out.psms,
             P_FRAGGER.out.psms)
    prot2intersect = P_COMET.out.prot2intersect
        .mix(P_IPY.out.prot2intersect,
             P_MSGF.out.prot2intersect,
             P_FRAGGER.out.prot2intersect)
    psm2combinedPEP = P_COMET.out.psm2combinedPEP
        .mix(P_IPY.out.psm2combinedPEP,
             P_MSGF.out.psm2combinedPEP,
             P_FRAGGER.out.psm2combinedPEP)
}
