process SMSNET {
    publishDir "$outdir", mode: 'symlink'
    conda "/home/shannc/anaconda3/envs/smsnet-shann"

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
    smsnet.py --model_dir $model_dir \
    --inference_input_file $mgfs \
    """
    //
}
