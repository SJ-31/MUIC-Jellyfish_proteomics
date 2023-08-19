include { COMET } from  '../modules/comet'
include { MSFRAGGER } from '../modules/msfragger'
include { PERCOLATOR as P_COMET } from '../modules/percolator'
include { PERCOLATOR as P_IPY} from '../modules/percolator'
include { PERCOLATOR as P_FRAGGER } from '../modules/percolator'
include { IDENTIPY } from '../modules/identipy'
include { MAKE_BK_DB } from '../modules/make_bk_db'

workflow bk_decoys {
    take:
    percolator_out
    seq_mapping
    header_mapping
    mzML_ch

    main:
    MAKE_BK_DB(percolator_out, header_mapping.first(), seq_mapping.first(),
               "$params.results/2-Second_pass/BK_databases")
    .flatten().filter( ~/.*\.fasta/ )
    .branch {
        comet: it.baseName =~ /comet*/
        identipy: it.baseName =~ /identipy*/
        msfragger: it.baseName =~ /msfragger*/
        }.set { bk_db }
    COMET(mzML_ch.collect(), "$params.results/2-Second_pass/Comet",
          bk_db.comet)
    IDENTIPY(mzML_ch.collect(), "$params.results/2-Second_pass/Identipy",
             bk_db.identipy)
    MSFRAGGER(mzML_ch.collect(), "$params.config/MSFragger_params.params",
    "$params.results/2-Second_pass/MsFragger",
              bk_db.msfragger)
    P_COMET(COMET.out.percolator,
            "$params.results/2-Second_pass/Percolator", bk_db.comet)
    P_IPY(IDENTIPY.out.percolator,
            "$params.results/2-Second_pass/Percolator", bk_db.identipy)
    P_FRAGGER(MSFRAGGER.out.percolator,
            "$params.results/2-Second_pass/Percolator", bk_db.msfragger)

    emit:
    prot2intersect = P_COMET.out.prot2intersect.mix(P_IPY.out.prot2intersect,
                                   P_FRAGGER.out.prot2intersect)
    psm2combinedPEP = P_COMET.out.psm2combinedPEP.mix(P_IPY.out.psm2combinedPEP,
                                   P_FRAGGER.out.psm2combinedPEP)
    prot2combinedPEP = P_COMET.out.prot2combinedPEP.mix(P_IPY.out.prot2combinedPEP,
                                   P_FRAGGER.out.prot2combinedPEP)
}
