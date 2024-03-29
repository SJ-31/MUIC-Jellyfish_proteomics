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
    Rscript $params.bin/R/sort_open_searches.r \
        -p . \
        -o temp.tsv \
        -r $params.bin/R \
        -m $seq_header_mappings

    Rscript $params.bin/R/unify_groups.r \
        -i temp.tsv \
        -o grouped_open_searches.tsv \
        -s open \
        -p O
    """
    //
}
