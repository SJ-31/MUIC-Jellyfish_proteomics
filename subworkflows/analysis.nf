analysis_logs = "$params.results/Analysis/Logs"

include { MAKE_ORGDB } from '../modules/make_org_db' addParams(logdir: analysis_logs)
include { ONTOLOGIZER } from '../modules/ontologizer' addParams(logdir: analysis_logs)
include { EMBEDDINGS } from '../modules/embeddings' addParams(logdir: analysis_logs)
include { EMBEDDINGS as GO_EMBEDDINGS } from '../modules/embeddings' addParams(logdir: analysis_logs)
include { COMPARISON_EMBEDDINGS } from '../modules/comparison_embeddings' addParams(logdir: analysis_logs)
include { GET_COG } from '../modules/get_cog' addParams(logdir: analysis_logs)
include { AGGREGATE } from '../modules/aggregate' addParams(logdir: analysis_logs, sem_distances: "/mnt/data/shannc/nf/results/C_indra/Analysis/sem_matrices.hdf5")
include { VIEW_ALIGNMENTS } from '../modules/view_alignments' addParams(logdir: analysis_logs)
include { GET_REPRESENTATIVES } from '../modules/get_representatives' addParams(logdir: analysis_logs)
include { DR } from '../modules/dimensionality_reduction' addParams(logdir: analysis_logs)
include { DR as DR_COMPARE } from '../modules/dimensionality_reduction' addParams(logdir: analysis_logs)

workflow 'analyze' {
    take:
    combined_tsv
    alignments
    peptide_map
    outdir

    main:
    MAKE_ORGDB(combined_tsv, outdir)
    GET_COG(combined_tsv, outdir)
    // VIEW_ALIGNMENTS(combined_tsv, alignments, peptide_map, "$outdir/Alignments")
    ONTOLOGIZER(combined_tsv, "$outdir/Ontologizer") // Overrepresentation analysis
    GO_EMBEDDINGS(combined_tsv, "$outdir/GO_term_embeddings", "a2v_go")
    model = "prottrans" // esm | prottrans | a2v (for GO embeddings)
    EMBEDDINGS(combined_tsv, "$outdir/Embeddings_${model}", model)
    COMPARISON_EMBEDDINGS(EMBEDDINGS.out.embd,
                          GET_COG.out.tsv,
                          params.comparison_prottrans,
                          params.uniprot_reviewed,
                          "protein",
                          "$outdir/Comparison_embeddings") // Model is esm or prottrans
    DR(EMBEDDINGS.out.embd, EMBEDDINGS.out.dist,
       GET_COG.out.tsv, "$params.config_dir/NO_FILE",
       "tsne", "", "$outdir/TSNE")
    DR_COMPARE(COMPARISON_EMBEDDINGS.out.embd,
               COMPARISON_EMBEDDINGS.out.dist,
               GET_COG.out.tsv, params.uniprot_reviewed,
               "tsne", "C", "$outdir/TSNE_comparison")
    // GET_REPRESENTATIVES(GO_PARENTS.out.categorized,
    //                     EMBEDDINGS.out.embd,
    //                     EMBEDDINGS.out.dist, outdir)
    // AGGREGATE(GO_PARENTS.out.categorized,
    //     EMBEDDINGS.out.embd,
    //     EMBEDDINGS.out.dist,
    //     "$outdir/Aggregated")
}
