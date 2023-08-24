process UNMATCHED_PSMS {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator)
    val(outdir)
    //

    output:
    path("unmatched_peptides.fasta")
    //

    script:
    """
    unmatched_peptides.py unmatched_peptides.fasta 0.05
    """
    //
}


process UNMATCHED_MSMS {
    publishDir "$outdir", mode: "copy"

    input:
    path(comet)
    path(identipy)
    path(msfragger)
    path(maxquant_all_pin)
    path(metamorpheus_all_psms)
    path(tide_target_search)
    path(msms_mapping)
    path(mzmls)
    val(outdir)
    //

    output:
    path("unmatched_*.mzML")
    //

    script:
    """
    Rscript $params.bin/unmatched_msms.r \
        --metamorpheus $metamorpheus_all_psms \
        --maxquant $maxquant_all_pin \
        --comet $comet \
        --identipy $identipy \
        --tide $tide_target_search \
        --msfragger $msfragger \
        -m $msms_mapping \
        --mzML_path .
    """
    //
}
