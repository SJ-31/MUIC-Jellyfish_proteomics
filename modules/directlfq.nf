process DIRECTLFQ {
    publishDir "$outdir", mode: "copy"

    input:
    path(comet)
    path(identipy)
    path(msfragger)
    path(maxquant_pin_file)
    path(metamorpheus_AllPSMs)
    path(tide_target_search)
    path(msms_mapping)
    val(outdir)
    //

    output:
    path("*")
    //

    script:
    """
    Rscript $params.bin/directlfq_format.r \
        -o directlfq.tsv \
        --metamorpheus $metamorpheus_AllPSMs \
        --maxquant $maxquant_pin_file \
        --comet $comet \
        --identipy $identipy \
        --tide $tide_target_search \
        --msfragger $msfragger \
        -m $msms_mapping \
        -t 0.05
    """
    //
    // directlfq lfq -i directlfq.tsv
}
