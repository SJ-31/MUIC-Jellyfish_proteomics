process METAMORPHEUS {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "${params.pref}_results*"
    memory "10 GB"
    stageInMode "copy"
    // debug true
    conda "/home/shannc/anaconda3/envs/metamorpheus"

    input:
    path(mzmls)
    val(outdir)
    val(mode)
    val(config)
    val(database)
    //

    output:
    path("metamorpheus${mode}*")
    path("metamorpheus${mode}_AllPSMs.psmtsv"), emit: psms
    tuple val("metamorpheus${mode}"), path("metamorpheus${mode}_AllPSMs_FormattedForPercolator.tab"), emit: percolator
    //

    shell:
    template 'metamorpheus.sh'
}
