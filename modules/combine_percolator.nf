process COMBINE_PERCOLATOR {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(percolator_protein)
    path(combined_results)
    path(seq_header_map)
    path(unmatched_peptides)
    val(outdir)
    //

    output:
    path("percolator_all.tsv")
    path("engine_stats.tsv")
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
        -s $seq_header_map \
        -u $unmatched_peptides \
        -r $params.bin/R

    cp .command.out get_percolator.log
    """
    //
}
