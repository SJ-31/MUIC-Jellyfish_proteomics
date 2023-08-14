include { COMET } from  '../modules/comet'
include { MSFRAGGER } from '../modules/msfragger'
include { PERCOLATOR } from '../modules/percolator'
include { IDENTIPY } from '../modules/identipy'
include { MAKE_BK_DB } from '../modules/make_bk_db'

workflow bk_decoys {
    take:
    percolator_out
    database
    mzXML_ch
    mzML_ch

    main:
    MAKE_BK_DB(percolator_out, params.mappings, database).flatten().filter( ~/.*\.fasta/ )
        .branch {
            comet: it[0] =~ /comet/
            identipy: it[0] =~ /identipy/
            msfragger: it[0] =~ /msfragger/
            }.set { bk_db }
    COMET(mzXML_ch.collect(), "$params.results/Second_pass/Comet",
          bk_db.comet.)
    IDENTIPY(mzML_ch.collect(), "$params.results/Second_pass/Identipy",
             bk_db.identipy)
    MSFRAGGER(mzML_ch.collect(), "$params.config/MSFragger_params.params",
    "$params.results/Second_pass/MsFragger",
              bk_db.msfragger)
}
