process SORT_OPEN {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_open_searches)
    path(seq_header_mappings)
    val(outdir)
    //

    output:
    path("open_searches.tsv")
    //

    script:
    """
    Rscript $params.bin/R/sort_open_searches.r \
        -p . \
        -o open_searches.tsv \
        -r $params.bin/R \
        -m $seq_header_mappings
    """
    //
}
