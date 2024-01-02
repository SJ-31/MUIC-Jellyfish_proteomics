process WRITE_QUANT {
    publishDir "$outdir", mode: "copy"

    input:
    path(directlfq)
    path(flashlfq)
    val(outdir)
    //

    output:
    path("sorted_directlfq.tsv"), emit: dlfq
    path("sorted_flashlfq.tsv"), emit: flfq
    //

    script:
    """
    write_quant.py \
        --dlfq $directlfq \
        --flfq $flashlfq \
        --dlfq_sorted sorted_directlfq.tsv \
        --flfq_sorted sorted_flashlfq.tsv
    """
    //
}
