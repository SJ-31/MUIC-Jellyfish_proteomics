process FALCON {
    publishDir "$outdir", mode: "copy"
    conda '/home/shannc/anaconda3/envs/falcon'


    input:
    path(mzMLs)
    //

    output:
    //

    shell:
    '''
    falcon ./*.mzML  \
        --export_representatives \
        falcon_cluster.csv

    find $(pwd) -name "*.mgf" > convert.txt
    msconvert -f convert.txt -o . --mzML
    find . \\( -name "*.mzML" -o -name "*.mgf" \\) \
        -type f > temp_list.txt
    cat temp_list.txt | sed 's;./;!{outdir}/;' > !{params.pref}.falcon.txt
    '''
    //
}
