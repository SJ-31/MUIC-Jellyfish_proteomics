process SORT_OPEN {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_open_searches)
    path(seq_header_mappings)
    val(outdir)
    //

    output:
    path("grouped_open_searches.tsv")
    //

    script:
    """
    Rscript $params.bin/sort_open_searches.r \
        -i . \
        -o grouped_open_searches.tsv \
        -r $params.bin \
        -m $seq_header_mappings
    """
    //
}
