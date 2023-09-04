process IDENTIPY {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/identipy"
    publishDir "$params.logs", mode: "copy", pattern: "*.log"

    input:
    path(mzMLs)
    val(outdir)
    val(database)
    //
    output:
    path("*.pepxml")
    tuple val("identipy"), path("identipy_all_pins.temp"), emit: percolator
    path("*.log")
    //
    shell:
    template 'identipy.sh'
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -method {reverse,shuffle} for decoy generation
    // -prefix decoy prefix
}
