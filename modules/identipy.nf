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
    def check = file("${outdir}/${mzML[0].baseName}.pep.xml")
    if (check.exists()) {
        '''
        cp !{outdir}/*.pep.xml .
        '''
    } else {
        template 'identipy.sh'
    }
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -method {reverse,shuffle} for decoy generation
    // -prefix decoy prefix
}

process FORMAT_IDPY {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/identipy"
    memory "35 GB"

    input:
    path(pepxmls)
    val(outdir)
    //

    output:
    tuple val("identipy"), path("identipy_all_pins.temp"), emit: percolator
    //

    shell:
    output = "${outdir}/identipy_all_pins.temp"
    if (file(output).exists()) {
        '''
        cp !{output} .
        '''
    } else {
        template 'format_identipy.sh'
    }
    //
}
