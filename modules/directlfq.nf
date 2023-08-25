process DIRECTLFQ_FORMAT {
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
    path("directlfq.aq_reformat.tsv")
    //

    script:
    """
    Rscript $params.bin/directlfq_format.r \
        -o directlfq.aq_reformat.tsv \
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
}

process DIRECTLFQ {
    conda '/home/shannc/anaconda3/envs/directlfq'
    publishDir "$outdir", mode: "copy"

    input:
    path(aqreformat)
    val(outdir)
    //

    output:
    path("directlfq_prot.tsv")
    //

    script:
    """
    directlfq lfq -i $aqreformat
    mv directlfq.aq_reformat.tsv.protein_intensities.tsv directlfq_prot.tsv
    """
    //
}
