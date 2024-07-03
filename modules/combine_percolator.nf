process COMBINE_PERCOLATOR {
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(percolator_protein)
    path(percolator_protein_open_search)
    path(combined_results)
    path(directlfq_input) // For mapping more peptide sequences
    path(seq_header_map)
    path(unmatched_peptides)
    val(outdir)
    val(logdir)
    //

    output:
    path("percolator_all.tsv")
    path("seq-header_map_found.tsv"), emit: seq_map
    path("percolator_peptide_map.tsv"), emit: peptide_map
    path("*.log")
    //

    script:
    """
    Rscript $params.bin/R/get_percolator_all.r \
        -p . \
        -o . \
        -i $combined_results \
        -d $directlfq_input \
        -s $seq_header_map \
        -u $unmatched_peptides \
        -r $params.bin/R

    cp .command.out get_percolator.log
    """
    //
}
