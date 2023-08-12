process SEARCH_INTERSECT {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_protein)
    val(outdir)
    //

    output:
    path("intersected_searches.tsv")
    //

    script:
    """
    Rscript $params.bin/atleast2.r \
        -m $params.mappings \
        -o intersected_searches.tsv
    """
    //
}
