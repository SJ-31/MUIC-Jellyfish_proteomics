process GO_PARENTS {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_results)
    val(outdir)
    //

    output:
    path("*.log")
    path("all_go_info.tsv")
    path("go_parents.json")
    path("missing_counts.tsv")
    path(output), emit: categorized
    //

    script:
    output = combined_results.baseName.replaceFirst(/wcoverage/, "wcategory") + ".tsv"
    """
    Rscript $params.bin/R/go_info.r -i $combined_results \
        -r $params.bin/R \
        -o . \
        -m info

    Rscript $params.bin/R/go_info.r --input $combined_results \
        --r_source $params.bin/R \
        --go_info_path all_go_info.tsv \
        --python_source $params.bin \
        --outdir . \
        --go_path $params.go \
        --n_groups 17 \
        --predefined_groups $params.go_predefined \
        --mode get_parents

    Rscript $params.bin/R/go_info.r --input $combined_results \
        --r_source $params.bin/R \
        --go_info_path all_go_info.tsv \
        --parent_mapping go_parents.json \
        --predefined_groups $params.go_predefined \
        --python_source $params.bin \
        --outdir . \
        --mode write_results

    cp .command.out go_parents.log
    """
    //
}
