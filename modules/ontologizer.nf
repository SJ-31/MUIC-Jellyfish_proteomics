process ONTOLOGIZER {
    // For overrepresentation analysis
    publishDir "$outdir", mode: "copy", pattern: "*.txt"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("*")
    path("ontologizer-*"), emit: over
    path("*_GO_slims.tsv")
    path("*.png")
    //

    script:
    """
    Rscript $params.bin/R/ontologizer.r \
        -m prep \
        -i $combined_tsv \
        --r_source $params.bin/R \
        --go_path $params.go \
        --python_source $params.bin \
        --executable $params.ontologizer_jar

    cp .command.log ontologizer.log

    Rscript $params.bin/R/ontologizer.r \
        -m get_slims \
        --results_path . \
        --r_source $params.bin/R \
        --go_slim_path $params.go_slims \
        --go_path $params.go

    Rscript $params.bin/R/ontologizer.r \
        -m word_cloud \
        --results_path . \
        --r_source $params.bin/R \
        --go_slim_path $params.go_slims \
        --go_path $params.go \
        --go_tm_dir $params.go_texts
    """
    //
}
