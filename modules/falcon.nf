process FALCON {
    publishDir "$outdir", mode: "copy"
    conda '/home/shannc/anaconda3/envs/falcon'


    input:
    path(mzMLs)
    val(outdir)
    //

    output:
    path("*.falcon.txt")
    path("*{mgf,mzML}")
    //

    shell:
    '''
    mkdir files; mv !{mzMLs} files
    falcon files/*.mzML  falcon \
        --export_representatives

    find $(pwd) -name "*.mgf" > convert.txt
    msconvert -f convert.txt -o . --mzML
    make_manifest.py -f . -t !{projectDir}/!{params.manifest_file} \
        -o !{params.pref}.falcon.txt -p !{outdir}
    '''
    //
}
