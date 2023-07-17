process SMSNET {
    publishDir "$outdir", mode: 'symlink'
    conda "/mnt/data/sirasris/miniconda3/envs/smsnet"

    input:
    path(mgfs)
    val(outdir)
    val(model_dir)
    //
    output:
    path("*")
    //
    script:
    """
    smsnet_wrapper.sh $mgfs $model_dir $outdir
    """
    //
}
