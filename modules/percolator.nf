process PERCOLATOR {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "*.log"

    input:
    tuple val(engine), path(pin_file)
    val(outdir)
    val(database)
    //

    output:
    tuple val(engine), path("${engine}_percolator_proteins.tsv"), path("${engine}_percolator_decoy_proteins.tsv"), emit: prot
    path("${engine}*.tsv")
    path("${engine}_percolator_proteins.tsv"), emit: prot2intersect
    path("${engine}_psm2combined_PEP.tsv"), emit: psm2combinedPEP
    path("${engine}_prot2combined_PEP.tsv"), emit: prot2combinedPEP
    path("*.log")
    //

    script:
    """
    percolator_wrapper_combined.sh \
        -p $engine \
        -i $pin_file \
        -f $database > percolator_${engine}.log

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
