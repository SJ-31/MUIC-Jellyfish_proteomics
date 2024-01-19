include { MAKE_ORGDB } from '../modules/make_org_db.nf'
include { ONTOLOGIZER } from '../modules/ontologizer.nf'

workflow 'analyze' {
    take:
    combined_tsv
    outdir

    main:
    MAKE_ORGDB(combined_tsv, outdir)
    ONTOLOGIZER(combined_tsv, outdir) // Overrepresentation analysis
}
