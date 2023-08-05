process SEARCH_INTERSECT {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_protein)
    //

    output:
    path("intersected_searches.tsv")
    //

    script:
    """
    Rscript $params.bin/atleast2.r \
        -m $params.mapping \
        -o intersected_searches.tsv
    """
    //
}
