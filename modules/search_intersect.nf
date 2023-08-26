process SEARCH_INTERSECT {
    publishDir "$outdir", mode: "copy"
    tag "Intersecting protein identifications: $percolator_protein"

    input:
    path(percolator_protein)
    val(outdir)
    val(header_mappings)
    val(seq_mappings)
    //

    output:
    path("unified_groups.tsv")
    path("intersected_searches.tsv"), emit: unsorted
    //

    script:
    """
    Rscript $params.bin/atleast2.r \
        -m $header_mappings \
        -o temp.tsv

    Rscript $params.bin/protein_coverage.r \
        -m $seq_mappings \
        -o intersected_searches.tsv

    Rscript $params.bin/unify_groups.r \
        intersected_searches.tsv \
        unified_groups.tsv
    """
    //
}
