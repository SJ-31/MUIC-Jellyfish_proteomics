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
    //

    script:
    """
    Rscript $params.bin/R/ontologizer.r \
        -m prep \
        -i $combined_tsv

    java -jar $params.ontologizer_jar -g $params.go \
            -a protein_mappings.ids \
            -p universe.txt \
            -s id_with_open.txt \
            -m "Bonferroni-Holm"

    java -jar $params.ontologizer_jar -g $params.go \
            -a protein_mappings.ids \
            -p universe.txt \
            -s unknown_to_db.txt \
            -m "Bonferroni-Holm"
    mv table-unknown_to_db* ontologizer-unknown_to_db.txt
    mv table-id_with_open* ontologizer-id_with_open.txt
    cp .command.log ontologizer.log

    Rscript $params.bin/R/ontologizer.r \
        -m get_slims \
        --results_path . \
        --go_slim_path $params.go_slims \
        --go_path $params.go
    """
    //
}
