process IDENTIPY {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/identipy"

    input:
    path(mzMLs)
    val(outdir)
    //
    output:
    path("${params.pref}_identipy.tsv")
    path("identipy*.tsv")
    path("${params.pref}-identipy_percolator")
    //
    shell:
    """
    identipy_wrapper.sh -d $params.databaseWdecoy \
        -p $params.pref \
        -c $params.bin/identipy2pin.r

    percolator_wrapper.sh \
        -p ${params.pref} \
        -i pin \
        -f ${params.database} \
        -e identipy
    """
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -method {reverse,shuffle} for decoy generation
    // -prefix decoy prefix
}
