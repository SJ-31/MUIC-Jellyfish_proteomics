process COMBINED_DATABASE {
    publishDir "$outdir", mode: "copy"

    input:
    path(other_fasta)
    path(denovo_peptidesWdecoys)
    path(denovo_peptides)
    path(denovo_decoys)
    val(outdir)
    //

    output:
    path("decoysWnormal.fasta")
    path("all_normal.fasta")
    path("all_decoys.fasta")
    //

    script:
    """
    make_combined.db $other_fasta
    cat decoysWnormal.temp $denovo_peptidesWdecoys > decoysWnormal.fasta
    cat all_decoys.temp $denovo_decoys > all_decoys.fasta
    cat all_normal.temp $denovo_peptides > all_normal.fasta
    """
    //
}
