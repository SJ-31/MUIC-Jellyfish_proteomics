process PEPNET {
    publishDir "$outdir", mode: "copy"
    conda "/mnt/data/shannc/anaconda3/envs/PepNet"

    input:
    path(mgf)
    path(outdir)
    //

    output:
    //

    script:
    """

    """
    //
}
