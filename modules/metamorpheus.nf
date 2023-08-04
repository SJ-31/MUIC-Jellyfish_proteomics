process METAMORPHEUS {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs/Metamorpheus", mode: "copy", pattern: "${params.pref}_results*"
    debug true
    conda "/home/shannc/anaconda3/envs/metamorpheus"

    input:
    path(mzmls)
    val(outdir)
    //

    output:
    path("metamorpheus*")
    tuple val("metamorpheus"), path("metamorpheus_AllPSMs_FormattedForPercolator.tab"), emit: percolator
    //

    shell:
    '''
    metamorpheus -s !{mzmls} \
        -o . \
        -t !{params.config}/metamorpheus_params.toml \
        -d !{params.database}

    mv Task1SearchTask/[Aa]ll* .

    for i in [Aa]ll*
        do
        mv $i metamorpheus_${i}
        done
    '''
}
