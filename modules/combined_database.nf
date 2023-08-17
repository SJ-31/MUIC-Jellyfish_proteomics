process COMBINED_DATABASE {
    publishDir "$outdir", mode: "copy"
    publishDir "$projectDir", mode: "copy", pattern: "*.DB.txt"

    input:
    path(other_fasta)
    path(denovo_peptides)
    val(outdir)
    //

    output:
    path("decoysWnormal.fasta")
    path("all_normal.fasta")
    path("all_decoys.fasta")
    path("*tsv")
    path("${params.pref}.DB.txt")
    //

    shell:
    '''
    outdir=!{outdir}
    make_combined_db.sh !{other_fasta} !{denovo_peptides}

    fasta_table.py decoysWnormal.fasta decoysWnormal_mapping.tsv
    fasta_table.py all_decoys.fasta all_decoys_mapping.tsv
    fasta_table.py all_normal.fasta all_normal_mapping.tsv
    find . -name "*[fasta,tsv]" -type f > temp_list.txt
    cat temp_list.txt | sed 's;./;${outdir}/;' > !{params.pref}.DB.txt
    '''
    //
}
