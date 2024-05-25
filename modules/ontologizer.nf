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
    groups=("id_with_open" "unknown_to_db" "unmatched_only")
    Rscript $params.bin/R/ontologizer.r \
        -m prep \
        -i $combined_tsv

    for group in "\${groups[@]}"; do
        java -jar $params.ontologizer_jar -g $params.go \
                -a protein_mappings.ids \
                -p universe.txt \
                -s "\$group".txt \
                -m "Bonferroni-Holm"
    done

    mv table-unknown_to_db* ontologizer-unknown_to_db.txt
    mv table-unmatched_only* ontologizer-unmatched_only.txt
    mv table-id_with_open* ontologizer-id_with_open.txt
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
