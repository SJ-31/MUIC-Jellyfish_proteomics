analysis_logs = "$params.results/Analysis/Logs"

include { MAKE_ORGDB } from '../modules/make_org_db' addParams(logdir: analysis_logs)
include { ONTOLOGIZER } from '../modules/ontologizer' addParams(logdir: analysis_logs)
include { EMBEDDINGS } from '../modules/embeddings' addParams(logdir: analysis_logs)
include { COMPARISON_EMBEDDINGS } from '../modules/comparison_embeddings' addParams(logdir: analysis_logs)

workflow 'analyze' {
    take:
    combined_tsv
    outdir

    main:
    MAKE_ORGDB(combined_tsv, outdir)
    ONTOLOGIZER(combined_tsv, "$outdir/Ontologizer") // Overrepresentation analysis
    EMBEDDINGS(combined_tsv, "$outdir/Embeddings", "esm") // Model is esm or prottrans
    COMPARISON_EMBEDDINGS(EMBEDDINGS.out.embd,
                          combined_tsv,
                          params.comparison_esm,
                          params.uniprot_reviewed,
                          "$outdir/Comparison_embeddings") // Model is esm or prottrans
}
