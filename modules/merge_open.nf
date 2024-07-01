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
    path(unmatched_peptides_out), emit: unmatched_pep
    //

    shell:
    def check = file("${outdir}/database_hits.tsv")
    unmatched_peptides_out = "unmatched_peptides-merge_open.tsv"
    if (check.exists()) {
        '''
        cp !{outdir}/unknown_hits.tsv .
        cp !{outdir}/unknown.fasta .
        cp !{outdir}/database_hits.tsv .
        cp !{outdir}/!{unmatched_peptides_out} .
        '''
    } else {
        '''
        merge_open.py \
            -i !{search_intersections} \
            --unmatched_peptides_in !{unmatched_peptides} \
            --unmatched_peptides_out !{unmatched_peptides_out} \
            --database_output database_hits.tsv \
            --unknown_output unknown_hits.tsv \
            --unknown_output_fasta temp1.fasta \
            -s !{open_searches}

        cat temp1.fasta unmatched_peptides.fasta > temp2.fasta
        cd-hit -i temp2.fasta -o unknown.fasta -c 1.0
        '''
    }
    // cd-hit at 1. identity will cluster sequences that are complete subsets
    // of one another
    //
}

// awk 'BEGIN {FS="\t";OFS="\t"} {if (FNR == 1) {print $0; next} else if($1 ~ /[TD]/) {print $0}}' \
