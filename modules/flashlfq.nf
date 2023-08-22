process FLASHFLQ {
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
    mkdir mzml; mv *.mzml mzml
    Rscript $params.bin/flashlfq_format.r \
        -o flashlfq.tsv \
        --metamorpheus $metamorpheus \
        --maxquant $maxquant \
        --comet $comet \
        --identipy $identipy \
        --tide $tide \
        --msfragger $msfragger \
        -m $msms_mapping

    $params.dotnet6 $params.flashlfq \
        --idt flashlfq.tsv \
        --out . \
        --rep mzml \
        --ppm 5 \
    """
    //
}
