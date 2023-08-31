process MERGE_QUANT {
    publishDir "$outdir", mode: "copy"

    input:
    path(directlfq)
    path(search_intersections)
    path(seq_mapping)
    val(outdir)
    //

    output:
    path("unified_groups_with_quant.tsv")
    //

    script:
    """
    merge_quantifications.py -d $directlfq \
        -i $search_intersections \
        -m $seq_mapping \
        -p 0.05 \
        -o unified_groups_with_quant.tsv
    """
    //
}
