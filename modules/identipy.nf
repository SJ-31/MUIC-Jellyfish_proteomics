process IDENTIPY {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/identipy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(mzML)
    val(outdir)
    val(logdir)
    val(database)
    //
    output:
    path("*.pep.xml"), emit: pepxml
    //
    shell:
    template 'identipy.sh'
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -method {reverse,shuffle} for decoy generation
    // -prefix decoy prefix
}

process FORMAT_IDPY {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/identipy"

    input:
    path(pepxmls)
    val(outdir)
    //

    output:
    tuple val("identipy"), path("identipy_all_pins.temp"), emit: percolator
    //

    shell:
    template 'format_identipy.sh'
    //
}
