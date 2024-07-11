process SMSNET {
    publishDir "$outdir", mode: 'copy'
    memory "10 GB"
    conda "/mnt/data/sirasris/miniconda3/envs/smsnet"

    input:
    path(mgf)
    val(outdir)
    //
    output:
    path(output)
    //
    shell:
    output = "${mgf.baseName}_smsnet.txt"
    check = file("${outdir}/${output}")
    if (check.exists()) {
        '''
        cp "!{outdir}/!{output}" .
        '''
    } else {
    '''
    export CUDA_VISIBLE_DEVICES=0
    python !{params.smsnet_exe} \
        --model_dir !{params.smsnetmodel} \
        --inference_input_file !{mgf}
    mv ./_output/!{mgf.baseName} ./!{mgf.baseName}_smsnet.txt
    '''
    }
   //
}

process EXTRACT_SMSNET {
    publishDir "$outdir", mode: "copy"

    input:
    path(smsnet_output)
    val(outdir)
    //

    output:
    path("smsnet_normal.fasta")
    //

    shell:
    '''
    cat !{smsnet_output} > smsnet_temp.txt
    filter_smsnet.sh
    '''
    //
}
