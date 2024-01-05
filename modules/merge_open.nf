process MERGE_OPEN {
    publishDir "$outdir", mode: "copy"

    input:
    path(search_intersections)
    path(open_searches) // Open searches need to be included as well, since they will map to de novo and transcriptome proteins
    path(unmatched_peptides)
    val(outdir)
    //

    output:
    path("unknown_hits.tsv"), emit: unknown_tsv
    path("unknown.fasta"), emit: unknown_fasta
    path("database_hits.tsv"), emit: database_tsv
    path("unmatched_peptides.tsv"), emit: unmatched_pep
    //

    shell:
    '''
    merge_open.py \
        -i !{search_intersections} \
        --unmatched_peptides !{unmatched_peptides} \
        --database_output database_hits.tsv \
        --unknown_output unknown_hits.tsv \
        --unknown_output_fasta temp1.fasta \
        -s !{open_searches}

    cat temp1.fasta unmatched_peptides.fasta > temp2.fasta
    cd-hit -i temp2.fasta -o unknown.fasta -c 1.0
    '''
    // cd-hit at 1.0 identity will cluster sequences that are complete subsets
    // of another
    //
}

// awk 'BEGIN {FS="\t";OFS="\t"} {if (FNR == 1) {print $0; next} else if($1 ~ /[TD]/) {print $0}}' \
