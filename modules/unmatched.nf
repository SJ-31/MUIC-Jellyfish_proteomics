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
    tuple val(engine), path(scan_file), path(percolator_file)
    path(msms_mapping)
    path(mzmls)
    val(outdir)
    //

    output:
    path("filtered_*.mzML")
    //

    script:
    """
    Rscript $params.bin/unmatched_msms.r \
        --engine $engine \
        --scan_file $scan_file \
        --percolator_file $percolator_file \
        --pep_thresh 1 \
        -m $msms_mapping \
        --mzML_path .
    """
    //
}
