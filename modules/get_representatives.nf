process GET_REPRESENTATIVES {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_tsv)
    path(embeddings)
    path(distances)
    val(outdir)
    //

    output:
    path("${params.pref}_all_representatives.tsv")
    path("*.log")
    //

    script:
    """
    clustering.py \
        -m $params.mmseqs \
        -i $combined_tsv \
        -o "${params.pref}_all_representatives_mmseqs.tsv"

    Rscript $params.bin/R/analysis \
        -i $combined_tsv \
        -o "${params.pref}_all_representatives.tsv"

    Rscript $params.bin/R/get_clusters.r \
            --r_source $params.bin/R \
            --combined_results "${params.pref}_all_representatives.tsv" \
            --sample_name $params.pref \
            --ontologizer_path $params.ontologizer_jar \
            --python_source $params.bin \
            --go_path $params.go \
            --embeddings $embeddings \
            --outdir . \
            --distances $distances

    cp .command.out get_representatives.log
    """
    //
}
