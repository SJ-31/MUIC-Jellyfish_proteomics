process DR {
    publishDir "$outdir", mode: "copy"
    publishDir "params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(embeddings)
    path(distances)
    path(combined_results)
    val(comparison_file)
    val(technique)
    val(compare)
    val(outdir)
    //

    output:
    path("*")
    //

    script:
    if (compare == "C") {
        """
        Rscript "$params.bin/R/analysis/visualize_dr.r" \
            -f . \
            -s $params.pref \
            -e $embeddings \
            -t $technique \
            -d $distances \
            -r "$params.bin/R" \
            -p $params.bin \
            -c $combined_results \
            -u $comparison_file \
            --compare
        """
    } else {
        """
        Rscript "$params.bin/R/analysis/visualize_dr.r" \
            -f . \
            -s $params.pref \
            -e $embeddings \
            -d $distances \
            -t $technique \
            -c $combined_results \
            -r "$params.bin/R" \
            -p $params.bin
        """
    }
    //
}
