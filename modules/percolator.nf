process PERCOLATOR {
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"
    memory "10 GB"


    input:
    tuple val(engine), path(pin_file)
    val(outdir)
    val(logdir)
    val(database)
    //

    output:
    tuple val(engine), path("${engine}_percolator_proteins.tsv"), path("${engine}_percolator_decoy_proteins.tsv"), emit: prot
    path("${engine}*.tsv")
    path("${engine}_percolator_psms.tsv"), emit: psms
    path("${engine}_percolator_proteins.tsv"), emit: prot2intersect
    path("${engine}_psm2combined_PEP.tsv"), emit: psm2combinedPEP
    //

    script:
    def check = file("${outdir}/${engine}_percolator_proteins.tsv")
    if (check.exists()) {
        """
        mv -Z ${outdir}/${engine}_percolator* .
        mv -Z ${outdir}/${engine}_psm2combined_PEP.tsv .
        """
    } else {
        """
        percolator_wrapper_combined.sh \
            -p $engine \
            -i $pin_file \
            -f $database

        Rscript $params.bin/R/2combinedPEP.r \
            -e $engine \
            -r $params.bin/R \
            -m ${engine}_percolator_psms.tsv \
            -d ${engine}_percolator_decoy_psms.tsv \
            -o ${engine}_psm2combined_PEP.tsv
        """
    }
    //
}

