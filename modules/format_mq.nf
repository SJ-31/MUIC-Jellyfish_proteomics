process FORMAT_MQ {
    publishDir "$outdir", mode: "copy"

    input:
    path(maxquant_msms)
    val(outdir)
    //

    output:
    tuple val("maxquant"), path("maxquant_all_pins.temp")
    //

    shell:
    template 'format_mq.sh'
    //
}
