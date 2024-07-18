process FALCON {
    publishDir "$outdir", mode: "copy"
    conda '/home/shannc/anaconda3/envs/falcon'


    input:
    path(mzMLs)
    val(outdir)
    //

    output:
    path("*.falcon.tsv")
    path("*{mgf,mzML}")
    //

    shell:
    def check = file("${outdir}/${params.pref}.falcon.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/* .
        '''
    } else {
    '''
    mkdir files; mv !{mzMLs} files
    falcon files/*.mzML  falcon \
        --export_representatives

    find $(pwd) -name "*.mgf" > convert.txt
    msconvert -f convert.txt -o . --mzML
    make_manifest.py -f . -t !{projectDir}/!{params.manifest_file} \
        -o !{params.pref}.falcon.tsv -p !{outdir}
    '''
    }

    //
}
