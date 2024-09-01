process ONTOLOGIZER {
    // For overrepresentation analysis
    publishDir "$outdir", mode: "copy", pattern: "{*.txt,*.png,*.tsv}"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("*")
    path("ontologizer-*.txt"), emit: over
    path("*_GO_slims.tsv")
    //

    script:
    """
    Rscript $params.bin/R/ontologizer.r \
        -m prep \
        -i $combined_tsv \
        --r_source $params.bin/R \
        --go_path $params.go \
        -w $projectDir \
        --python_source $params.bin \
        --executable $params.ontologizer_jar

    cp .command.log ontologizer.log

    Rscript $params.bin/R/ontologizer.r \
        -m get_slims \
        --results_path . \
        -w $projectDir \
        --r_source $params.bin/R \
        --go_slim_path $params.go_slims \
        --go_path $params.go
    """
    //
}
