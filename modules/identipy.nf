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
        -db $params.database \
        -of tsv \
        -at \
        -ad \
        -out .
    """
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -ad Add decoy
    // -method {reverse,shuffle} for decoy generation
    //
}
