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
    path("metamorpheus${mode}*"), emit: all
    path("metamorpheus${mode}_AllPSMs.psmtsv"), emit: psms, optional: true
    tuple val("metamorpheus${mode}"), path("metamorpheus${mode}_AllPSMs_FormattedForPercolator.tab"), emit: percolator, optional: true
    //

    shell:
    def check = file("${outdir}/{metamorpheus${mode}_AllPSMs*,*xml}")
    if (check && config =~ /gptmd/) { // For the "GPTMD" option, which prepares
        // an xml file
        '''
        mv -Z !{outdir}/*xml .
        '''
    } else if (check) {
        '''
        mv -Z !{outdir}/*{tsv,tab,txt} .
        '''
    } else {
        template 'metamorpheus.sh'
    }
}
