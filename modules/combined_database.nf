process COMBINED_DATABASE {
    publishDir "$outdir", mode: "copy"
    publishDir "$text_loc", mode: "copy", pattern: "*.DB.txt"

    input:
    path(other_fasta)
    path(denovo_peptides)
    val(outdir)
    val(text_loc)
    //

    output:
    path("decoysWnormal.fasta")
    path("all_normal.fasta")
    path("all_decoys.fasta")
    path("blast_db")
    path("*tsv")
    path("${params.pref}.DB.txt"), emit: listing
    //

    shell:
    '''
    outdir=!{outdir}
    make_combined_db.sh !{other_fasta} !{denovo_peptides}

    fasta_table.py decoysWnormal.fasta decoysWnormal_mapping.tsv
    fasta_table.py all_decoys.fasta all_decoys_mapping.tsv
    fasta_table.py all_normal.fasta all_normal_mapping.tsv
    find . \\( -name "*fasta" -o -name "*tsv" \\) -type f > temp_list.txt
    cat temp_list.txt | sed 's;./;!{outdir}/;' > !{params.pref}.DB.txt

    mv downloaded.fasta temp; cat temp | seqkit rmdup > downloaded.fasta
    makeblastdb -in downloaded.fasta \
        -parse_seqids \
        -blastdb_version 5 \
        -title Cnidarian_proteins \
        -dbtype prot \
        -out downloaded
    mkdir blast_db; mv downloaded* blast_db
    '''
    //
}
