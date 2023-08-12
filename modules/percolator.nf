process PERCOLATOR {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(engine), path(pin_file)
    val(outdir)
    //

    output:
    path("${engine}*.tsv")
    path("${engine}_percolator_proteins.tsv"), emit: prot2intersect
    path("${engine}_psm2combined_PEP.tsv"), emit: psm2combinedPEP
    path("${engine}_prot2combined_PEP.tsv"), emit: prot2combinedPEP
    //

    script:
    """
    percolator_wrapper_combined.sh \
        -p $engine \
        -i $pin_file \
        -f $params.database

    Rscript $params.bin/to_combined_PEP.r \
        -m ${engine}_percolator_psms.tsv \
        -d ${engine}_percolator_decoy_psms.tsv \
        -o ${engine}_psm2combined_PEP.tsv

    Rscript $params.bin/to_combined_PEP.r \
        -m ${engine}_percolator_proteins.tsv \
        -d ${engine}_percolator_decoy_proteins.tsv \
        --protein_matches \
        -o ${engine}_prot2combined_PEP.tsv
    """
    //
}
