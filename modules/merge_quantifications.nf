process MERGE_QUANT {
    publishDir "$outdir", mode: "copy"

    input:
    path(directlfq)
    path(search_intersections)
    val(outdir)
    //

    output:
    path("unknown_hits.tsv"), emit: unknown_tsv
    path("unknown.fasta"), emit: unknown_fasta
    path("database_hits.tsv"), emit: database
    //

    shell:
    '''
    merge_quantifications.py -d !{directlfq} \
        -i !{search_intersections} \
        -p !{params.pep_thresh} \
        -o database_hits.tsv

    awk -F "\t" '{if (FNR == 1) {print $0; next} else {$1 ~ /P/} {print $0}}' \
        temp.tsv > database_hits.tsv
    awk -F "\t" '{if (FNR == 1) {print $0; next} else {$1 ~ /O/} {print $0}}' \
        temp.tsv > unknown_hits.tsv
    awk -F "\t" '{printf ">%s\n%s\n",$1,$7}' unknown.fasta
    '''
    //
}
