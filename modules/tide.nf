process TIDE {
    publishDir "$outdir", mode: "copy", pattern: "tide-search*"
    publishDir "$percolatordir", mode: "copy", pattern: "*percolator*"
    publishDir "$logdir", mode: "copy", pattern: "*.log.txt"
    memory "35 GB"

    input:
    path(mgf)
    val(outdir)
    val(logdir)
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
    def check = file("${outdir}/tide-search.target.txt")
    if (check.exists()) {
        """
        cp ${outdir}/tide-search* .
        cp ${percolatordir}/tide_percolator* .
        """
    } else {
        template 'tide.sh'
    }
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
    //

    script:
    """
    Rscript $params.bin/2combinedPEP.r \
        -e tide \
        -m tide_percolator_psms.tsv \
        -d tide_percolator_decoy_psms.tsv \
        -o tide_psm2combined_PEP.tsv
    """
    //
}
