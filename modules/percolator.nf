process PERCOLATOR {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(engine), path(pin_file)
    //

    output:
    path("${engine}*.tsv")
    //

    script:
    """
    percolator_wrapper_combined.sh \
        -p $engine \
        -i $pin_file \
        -f $params.database

    Rscript $params.bin/to_combined_PEP.r \
        -m comet_percolator_psms.tsv \
        -d comet_percolator_decoy_psms.tsv \
        -o comet_psm2combined_PEP.tsv

    Rscript $params.bin/to_combined_PEP.r \
        -m comet_percolator_proteins.tsv \
        -d comet_percolator_decoy_proteins.tsv \
        --protein_matches \
        -o comet_prot2combined_PEP.tsv
    """
    //
}
