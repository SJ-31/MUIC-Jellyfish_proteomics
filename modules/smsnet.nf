process SMSNET {
    publishDir "$outdir", mode: 'copy'
    conda "/mnt/data/sirasris/miniconda3/envs/smsnet"

    input:
    path(mgfs)
    val(outdir)
    //
    output:
    path("${mgf.baseName}_smsnet.txt")
    //
    shell:
    '''
    export CUDA_VISIBLE_DEVICES=0
    python !{params.smsnet_exe} \
        --model_dir !{params.smsnetmodel} \
        --inference_input_file !{mgf} \
        --rescore
    mv ._output/!{mgf.baseName} ./!{mgf.baseName}_smsnet.txt
    '''
    //
}

process EXTRACT_SMSNET {
    publishDir "$outdir", mode: "copy"

    input:
    path(smsnet_output)
    val(outdir)
    //

    output:
    path(smsnet_normal.fasta)
    //

    shell:
    '''
    cat !{smsnet_output} > smsnet_temp.txt
    grep "<s>" -v smsnet_temp.txt | grep -v "<unk>" | grep -v "</s>" | \
        sed -e 's/ //g' -e 's/m//g' | \
        awk 'length >= 7 {printf ">SMSNet%s-DENOVO\n%s\n",FNR,$0}' | \
        > smsnet_normal.fasta
    '''
    //
}
