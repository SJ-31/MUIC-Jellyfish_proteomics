process MERGE_QUANT {
    publishDir "$outdir", mode: "copy"

    input:
    path(directlfq)
    path(search_intersections)
    path(seq_mapping)
    val(outdir)
    //

    output:
    path("final_protein_groups.tsv")
    //

    script:
    """
    merge_quantifications.py -d $directlfq \
        -i $search_intersections \
        -m $seq_mapping \
        -o final_protein_groups.tsv
    """
    //
}
