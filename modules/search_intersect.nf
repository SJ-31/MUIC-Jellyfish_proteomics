process SEARCH_INTERSECT {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_protein)
    val(outdir)
    val(header_mappings)
    //

    output:
    path("intersected_searches.tsv")
    //

    script:
    """
    Rscript $params.bin/atleast2.r \
        -m $header_mappings \
        -o intersected_searches.tsv
    """
    //
}
