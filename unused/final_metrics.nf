process FINAL_METRICS {
    // Calculate coverage, emPAI and
    publishDir "$outdir", mode: "copy"

    input:
    val(combined_table)
    val(outdir)
    //

    output:
    //
    // TODO: Unfinished, this needs to
    script:
    """
    Rscript $params.bin/R/protein_coverage.r \
        --intersected_searches $combined_table \
        -o ${params.pref}_final.tsv

    """
    //
}
