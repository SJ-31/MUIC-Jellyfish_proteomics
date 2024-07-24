process MERGE_OPEN {

    publishDir "$outdir", mode: "copy"

    input:
    path(search_intersections)
    path(open_searches) // Open searches need to be included as well, since they will map to de novo and transcriptome proteins
    val(outdir)
    //

    output:
    path("unknown_hits.tsv"), emit: unknown_tsv
    path("unknown_*.fasta"), emit: unknown_fasta
    path("database_hits.tsv"), emit: database_tsv
    path("query_map.tsv"), emit: query_tsv
    //

    shell:
    def check = file("${outdir}/database_hits.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/unknown_hits.tsv .
        cp !{outdir}/unknown_*.fasta .
        cp !{outdir}/database_hits.tsv .
        cp !{outdir}/query_map.tsv .
        '''
    } else {
        '''
        merge_open.py \
            -i !{search_intersections} \
            --database_output database_hits.tsv \
            --unknown_output unknown_hits.tsv \
            --unknown_output_fasta temp.fasta \
            -q query_map.tsv \
            -s !{open_searches}

        for i in temp*.fasta; do
            new_name=$(echo $i | sed s'/temp/unknown/')
            cd-hit -i $i -o $new_name -c 1.0
        done
        '''
    }
    // cd-hit at 1. identity will cluster sequences that are complete subsets
    // of one another
    //
}

// awk 'BEGIN {FS="\t";OFS="\t"} {if (FNR == 1) {print $0; next} else if($1 ~ /[TD]/) {print $0}}' \
