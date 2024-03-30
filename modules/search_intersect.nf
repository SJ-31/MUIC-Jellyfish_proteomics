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
    def check = file("${outdir}/unified_groups.tsv")
    if (check.exists()) {
        """
        mv -Z ${outdir}/unified_groups.tsv .
        mv -Z ${outdir}/intersected_searches.tsv .
        """
    } else {
        """
        Rscript $params.bin/R/atleast2.r \
            -m ${seq_header_mappings} \
            -p . \
            -r $params.bin/R \
            -o intersected_searches.tsv

        Rscript $params.bin/R/unify_groups.r \
            -i intersected_searches.tsv \
            -o unified_groups.tsv \
            -s standard \
            -p G
        """
    }
    // Leave the filtering for the combining module
}
