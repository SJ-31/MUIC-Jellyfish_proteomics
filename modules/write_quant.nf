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
    check = "sorted_directlfq.tsv"
    if (file("${outdir}/${check}").exists()) {
        """
        cp ${outdir}/sorted_directlfq.tsv .
        cp ${outdir}/sorted_flashlfq.tsv .
        """
    } else {
    """
    write_quant.py \
        --dlfq $directlfq \
        --flfq $flashlfq \
        --dlfq_sorted sorted_directlfq.tsv \
        --flfq_sorted sorted_flashlfq.tsv
    """
    }
    //
}
