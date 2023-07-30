process IDENTIPY {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/identipy"

    input:
    path(mzML)
    val(outdir)
    //
    output:
    path("${mzML.baseName}_identipy.tsv")
    //
    script:
    """
    identipy $mzML \
        -db $params.databaseWdecoy \
        -of tsv \
        -at \
        -prefix $params.decoy_prefix
        -out .
    mv ${mzML.baseName}.tsv ${mzML.baseName}_identipy.tsv
    Rscript $params.bin/identipy2pin.r ${mzML.baseName}_identipy.tsv \
        ${mzML.baseName}_identipy.pin \
        ${mzML.baseName}
    """
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -method {reverse,shuffle} for decoy generation
    // -prefix decoy prefix
}
