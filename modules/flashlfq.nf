process FLASHLFQ {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "*{txt,toml}"
    errorStrategy 'ignore'

    input:
    path(comet)
    path(identipy)
    path(msfragger)
    path(maxquant_pin_file)
    path(metamorpheus_AllPSMs)
    path(tide_target_search)
    path(msms_mapping)
    path(mzmls)
    val(outdir)
    //

    output:
    path("Quantified*")
    path("*.{txt,toml}")
    //

    script:
    """
    mkdir mzML; mv *.mzML mzML
    Rscript $params.bin/flashlfq_format.r \
        -o flashlfq.tsv \
        --metamorpheus $metamorpheus_AllPSMs \
        --maxquant $maxquant_pin_file \
        --comet $comet \
        --identipy $identipy \
        --tide $tide_target_search \
        --msfragger $msfragger \
        -m $msms_mapping \
        -t 0.05

    $params.dotnet6 $params.flashlfq \
        --idt flashlfq.tsv \
        --out . \
        --rep mzML \
        --ppm 5 \
    """
    //
}
