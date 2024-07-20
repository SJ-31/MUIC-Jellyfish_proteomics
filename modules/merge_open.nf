process MERGE_OPEN {

    publishDir "$outdir", mode: "copy"

    input:
    path(search_intersections)
    path(open_searches) // Open searches need to be included as well, since they will map to de novo and transcriptome proteins
    val(outdir)
    //

    output:
    path("unknown_hits.tsv"), emit: unknown_tsv
    path("unknown.fasta"), emit: unknown_fasta
    path("database_hits.tsv"), emit: database_tsv
    //

    shell:
    def check = file("${outdir}/database_hits.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/unknown_hits.tsv .
        cp !{outdir}/unknown.fasta .
        cp !{outdir}/database_hits.tsv .
        '''
    } else {
        '''
        merge_open.py \
            -i !{search_intersections} \
            --database_output database_hits.tsv \
            --unknown_output unknown_hits.tsv \
            --unknown_output_fasta temp1.fasta \
            -s !{open_searches}

        cd-hit -i temp1.fasta -o unknown.fasta -c 1.0
        '''
    }
    // cd-hit at 1. identity will cluster sequences that are complete subsets
    // of one another
    //
}

// awk 'BEGIN {FS="\t";OFS="\t"} {if (FNR == 1) {print $0; next} else if($1 ~ /[TD]/) {print $0}}' \
