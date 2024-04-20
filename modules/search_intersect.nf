process SEARCH_INTERSECT {
    publishDir "$outdir", mode: "copy"
    tag "Intersecting protein identifications: $percolator_protein"

    input:
    path(percolator_protein)
    val(outdir)
    val(seq_header_mappings)
    //

    output:
    path("intersected_searches.tsv"), emit: unsorted
    //

    script:
    def check = file("${outdir}/intersected_searches.tsv")
    if (check.exists()) {
        """
        cp ${outdir}/intersected_searches.tsv .
        """
    } else {
        """
        Rscript $params.bin/R/atleast2.r \
            -m ${seq_header_mappings} \
            -p . \
            -r $params.bin/R \
            -o intersected_searches.tsv \
            -s standard
        """
    }
    // Leave the filtering for the combining module
}
