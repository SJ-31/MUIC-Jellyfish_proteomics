process COMBINED_DATABASE {
    publishDir "$outdir", mode: "copy"
    publishDir "$mapdir", mode: "copy", pattern: "*tsv"

    input:
    path(other_fasta)
    path(denovo_peptidesWdecoys)
    path(denovo_peptides)
    path(denovo_decoys)
    val(outdir)
    val(mapdir)
    //

    output:
    path("decoysWnormal.fasta")
    path("all_normal.fasta")
    path("all_decoys.fasta")
    path("*tsv")
    //

    script:
    """
    make_combined.db $other_fasta
    cat decoysWnormal.temp $denovo_peptidesWdecoys | seqkit rmdup \
        > decoysWnormal.fasta
    cat all_decoys.temp $denovo_decoys | seqkit rmdup \
        > all_decoys.fasta
    cat all_normal.temp $denovo_peptides | seqkit rmdup \
        > all_normal.fasta
    fasta_table.py decoysWnormal.fasta decoysWnormal_mapping.tsv
    fasta_table.py all_decoys.fasta all_decoys_mapping.tsv
    fasta_table.py all_normal.fasta all_normal_mapping.tsv
    """
    //
}
