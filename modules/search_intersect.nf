process SEARCH_INTERSECT {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_protein)
    val(outdir)
    val(header_mappings)
    //

    output:
    path("unified_groups.tsv")
    path("intersected_searches.tsv"), emit: unsorted
    //

    script:
    """
    Rscript $params.bin/atleast2.r \
        -m $header_mappings \
        -o intersected_searches.tsv

    Rscript $params.bin/unify_groups.r \
        intersected_searches.tsv \
        unified_groups.tsv
    """
    //
}
