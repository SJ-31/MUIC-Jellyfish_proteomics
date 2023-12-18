process MERGE_QUANT {
    publishDir "$outdir", mode: "copy"

    input:
    path(directlfq)
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
    '''
    merge_quantifications.py -d !{directlfq} \
        -i !{search_intersections} \
        -s !{open_searches} \
        -p !{params.pep_thresh} \
        -q !{params.fdr} \
        -o temp.tsv

    awk 'BEGIN {FS="\t";OFS="\t"} {if (FNR == 1) {print $0; next} else if($1 ~ /P/) {print $0}}' \
        temp.tsv > database_hits.tsv
    awk 'BEGIN {FS="\t";OFS="\t"} {if (FNR == 1) {print $0; next} else if($1 ~ /O/) {print $0}}' \
        temp.tsv > unknown_hits.tsv
    awk -F "\t" '{if (FNR == 1) {next} else {printf ">%s\\n%s\\n",$1,$7}}' \
        unknown_hits.tsv > temp.fasta
    cd-hit -i temp.fasta -o unknown.fasta -c 1.0
    '''
    // cd-hit at 1.0 identity will cluster sequences that are complete subsets
    // of another
    //
}

    // awk 'BEGIN {FS="\t";OFS="\t"} {if (FNR == 1) {print $0; next} else if($1 ~ /[TD]/) {print $0}}' \
