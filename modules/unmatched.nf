process UNMATCHED_PSMS {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator)
    val(outdir)
    //

    output:
    path("unmatched_peptides.tsv"), emit: tsv
    // Obtain the tsv for sorting between high- and low-confidence unmatched
    // peptides

    script:
    """
    unmatched_peptides.py -i "." \
        -t unmatched_peptides.tsv \
        -p 1 \
        -q 0.05
    """
    //
}


process UNMATCHED_MSMS {
    publishDir "$outdir", mode: "copy"

    input:
    path(scan_mappings)
    path(psm2combinedPEP)
    path(mzmls)
    val(outdir)
    //

    output:
    path("filtered_*.mzML")
    //
    // Use the psm2combinedPEP csv file for this
    script:
    """
    mkdir psms; mv $psm2combinedPEP psms
    mkdir scans; mv $scan_mappings scans
    Rscript ${params.bin}/R/unmatched_msms.r \
        --psm_path psms  \
        --pep_thresh 1 \
        -r $params.bin/R \
        --scan_path scans \
        --mzML_path . > unmatched.txt
    """
    // The pep_thresh argument allows you to filter out msms spectra that gave rise to psms falling that fell below the given PEP threshold
}
