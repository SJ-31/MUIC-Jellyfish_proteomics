process GO_REPRESENTATIVES {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_results)
    val(outdir)
    //

    output:
    path("*.log")
    path("all_go_info.tsv")
    path("go_representatives.json")
    //

    script:
    """
    Rscript $params.bin/go_info.r -i $combined_results \
        -r $params.bin/R \
        -o .

    go_subset.py --go_path params.go \
        --sample_path $combined_results \
        --go_info_path all_go_info.tsv \
        -n 18 \
        --outdir .
    """
    //
}
