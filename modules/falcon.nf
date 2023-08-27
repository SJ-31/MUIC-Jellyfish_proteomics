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
    falcon !{mzMLs}  \
        --export_representatives \
        falcon_cluster.csv

    find $(pwd) -name "*.mgf" > convert.txt
    msconvert -f convert.txt -o . --mzML
    find . \\( -name "*.mzML" -o -name "*.mgf" \\) \
        -type f > temp_list.txt
    make_manifest.py -f . -t !{params.manifest_file} \
    -o !{params.pref}.falcon.txt
    '''
    //
}
