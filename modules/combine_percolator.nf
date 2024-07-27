process COMBINE_PERCOLATOR {
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(percolator_protein)
    path(percolator_protein_open_search)
    path(combined_results)
    path(directlfq_input) // For mapping more peptide sequences
    path(seq_header_map)
    val(outdir)
    val(logdir)
    //

    output:
    path("percolator_all.tsv")
    path("seq-header_map_found.tsv"), emit: seq_map
    path("percolator_peptide_map.tsv"), emit: peptide_map
    path("*.log"), optional: true
    //

    script:
    def check = file("${outdir}/percolator_all.tsv")
    if (check.exists()) {
        """
        cp $check .
        cp ${outdir}/seq-header_map_found.tsv .
        cp ${outdir}/percolator_peptide_map.tsv .
        """
    } else {
    """
    Rscript $params.bin/R/get_percolator_all.r \
        -p . \
        -o . \
        -i $combined_results \
        -d $directlfq_input \
        -s $seq_header_map \
        -r $params.bin/R \
        -y $params.bin

    cp .command.out get_percolator.log
    """
    }
    //
}
