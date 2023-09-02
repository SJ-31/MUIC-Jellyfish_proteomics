process TIDE {
    publishDir "$outdir", mode: "copy", pattern: "tide-search*"
    publishDir "$percolatordir", mode: "copy", pattern: "*percolator*"
    publishDir "$params.logs", mode: "copy", pattern: "*.log.txt"

    input:
    path(mgf)
    val(outdir)
    val(percolatordir)
    val(database)
    //

    output:
    path("tide-search.target.txt"), emit: target
    path("tide-search.decoy.txt"), emit: decoy
    path("*percolator*"), emit: percolator
    path("tide_percolator_proteins.tsv"), emit: perc_protein
    path("tide_percolator_psms.tsv"), emit: perc_psms
    path("*.log.txt")
    //

    script:
    template 'tide.sh'
    //
}


process TIDE_COMBINED_PEP {
    publishDir "$outdir", mode: "copy"

    input:
    path(percolator_files)
    val(outdir)
    //

    output:
    path("tide_psm2combined_PEP.tsv"), emit: psm2combinedPEP
    path("tide_prot2combined_PEP.tsv"), emit: prot2combinedPEP
    //

    script:
    """
    Rscript $params.bin/to_combined_PEP.r \
        -m tide_percolator_psms.tsv \
        -d tide_percolator_decoy_psms.tsv \
        -o tide_psm2combined_PEP.tsv

    Rscript $params.bin/to_combined_PEP.r \
        -m tide_percolator_proteins.tsv \
        -d tide_percolator_decoy_proteins.tsv \
        --protein_matches \
        -o tide_prot2combined_PEP.tsv
    """
    //
}
