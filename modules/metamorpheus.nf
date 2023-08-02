process METAMORPHEUS {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs/Metamorpheus", mode: "copy", pattern: "${params.pref}_results*"
    conda "/home/shannc/anaconda3/envs/metamorpheus"

    input:
    path(mzmls)
    val(outdir)
    //

    output:
    path("${params.pref}")
    path("${params.pref}-comet_percolator")
    //

    shell:
    '''
    metamorpheus -s !{mzmls} \
        -o . \
        -t !{params.morpheusPars} \
        -d !{params.databaseWdecoy} \
    mv Task1SearchTask/[Aa]ll* .
    percolator_wrapper.sh \
        -p !{params.pref} \
        -i AllPSMS_FormattedForPercolator.tab \
        -f !{params.database} \
        -e metamorpheus
    for i in [Aa]ll*
        do
        mv $i !{params.pref}_${i}
        done
    '''
}
