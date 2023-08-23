process DIRECTLFQ {
    publishDir "$outdir", mode: "copy"

    input:
    path(comet)
    path(identipy)
    path(msfragger)
    path(maxquant)
    path(metamorpheus)
    path(tide)
    path(msms_mapping)
    path(mzmls)
    val(outdir)
    //

    output:
    //

    script:
    """
    Rscript $params.bin/directlfq_format.r \
        -o flashlfq.tsv \
        --metamorpheus $metamorpheus \
        --maxquant $maxquant \
        --comet $comet \
        --identipy $identipy \
        --tide $tide \
        --msfragger $msfragger \
        -m $msms_mapping
    """
    //
}
