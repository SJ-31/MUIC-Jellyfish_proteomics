process SEARCH_INTERSECT {
    publishDir "$outdir", mode: "copy"
    tag "Intersecting protein identifications: $percolator_protein"

    input:
    path(percolator_protein)
    val(outdir)
    val(seq_header_mappings)
    //

    output:
    path("unified_groups.tsv"), emit: sorted
    path("intersected_searches.tsv"), emit: unsorted
    //

    script:
    """
    Rscript $params.bin/atleast2.r \
        -m ${seq_header_mappings} \
        -o intersected_searches.tsv \
        -p $params.pep_thresh \
        -f $params.fdr

    Rscript $params.bin/unify_groups.r \
        -i intersected_searches.tsv \
        -o unified_groups.tsv \
        -s standard \
        -p G
    """
    //
}
