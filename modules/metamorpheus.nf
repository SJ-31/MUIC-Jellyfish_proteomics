process METAMORPHEUS {
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "${params.pref}_results*"
    memory "10 GB"
    stageInMode "copy"
    // debug true
    conda "/home/shannc/anaconda3/envs/metamorpheus"

    input:
    path(mzmls)
    val(outdir)
    val(logdir)
    val(mode)
    val(config)
    val(database)
    //

    output:
    path("metamorpheus${mode}*")
    path("metamorpheus${mode}_AllPSMs.psmtsv"), emit: psms, optional: true
    tuple val("metamorpheus${mode}"), path("metamorpheus${mode}_AllPSMs_FormattedForPercolator.tab"), emit: percolator, optional: true
    //

    shell:
    template 'metamorpheus.sh'
}
