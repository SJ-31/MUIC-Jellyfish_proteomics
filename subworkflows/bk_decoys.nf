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
    database
    mzXML_ch
    mzML_ch

    main:
    MAKE_BK_DB(percolator_out, params.mappings, database,
               "$params.results/Second_pass/BK_databases")
    .flatten().filter( ~/.*\.fasta/ )
    .branch {
        comet: it.baseName =~ /comet*/
        identipy: it.baseName =~ /identipy*/
        msfragger: it.baseName =~ /msfragger*/
        }.set { bk_db }
    bk_db.identipy.view()
    COMET(mzXML_ch.collect(), "$params.results/Second_pass/Comet",
          bk_db.comet)
    IDENTIPY(mzML_ch.collect(), "$params.results/Second_pass/Identipy",
             bk_db.identipy)
    MSFRAGGER(mzML_ch.collect(), "$params.config/MSFragger_params.params",
    "$params.results/Second_pass/MsFragger",
              bk_db.msfragger)
    P_COMET(COMET.out.percolator,
            "$params.results/Second_pass/Percolator", bk_db.comet)
    P_IPY(IDENTIPY.out.percolator,
            "$params.results/Second_pass/Percolator", bk_db.identipy)
    P_FRAGGER(MSFRAGGER.out.percolator,
            "$params.results/Second_pass/Percolator", bk_db.msfragger)
}