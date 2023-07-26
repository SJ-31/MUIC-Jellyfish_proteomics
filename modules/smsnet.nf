process SMSNET {
    publishDir "$outdir", mode: 'copy'
    conda "/mnt/data/sirasris/miniconda3/envs/smsnet"

    input:
    path(mgfs)
    val(outdir)
    //
    output:
    path("*")
    //
    script:
    """
    smsnet_wrapper.sh $mgfs $params.smsnetmodel $outdir
    """
    //
}
