process AGGREGATE {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_results)
    path(embeddings)
    path(distances)
    val(outdir)
    //

    output:
    path("cluster-aggregated.tsv")
    path("cluster-aggregated_slims.tsv")
    path("Group-aggregated.tsv")
    path("Group-aggregated_slims.tsv")
    path("*.log")
    //

    script:
    """
    Rscript $params.bin/R/get_clusters.r \
            --r_source $params.bin/R \
            --combined_results $combined_results \
            --sample_name $params.pref \
            --ontologizer_path $params.ontologizer_jar \
            --python_source $params.bin \
            --go_path $params.go \
            --embeddings $embeddings \
            --outdir . \
            --distances $distances

    cp .command.out aggregate.log
    """
    //
}
