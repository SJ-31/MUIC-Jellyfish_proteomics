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


process FILTER_MSMS {
    publishDir "$outdir", mode: "copy"

    input:
    path(psm2combinedPEP)
    path(scan_mapping)
    path(mzmls)
    val(outdir)
    //

    output:
    path("filtered_*.mzML")
    //
    // Use the prot2combinedPEP csv file for this
    shell:
    '''
    Rscript !{params.bin}/unmatched_msms.r \
        --engine !{engine} \
        --scan_file all_scans.txt \
        --percolator_file $percolator_file \
        --pep_thresh 1 \
        -m $msms_mapping \
        --mzML_path . > unmatched.txt
    '''
    // The pep_thresh argument allows you to filter out msms spectra that gave rise to psms falling that fell below the given PEP threshold
}
