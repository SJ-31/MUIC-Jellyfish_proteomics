process XTANDEM {
    publishDir "$outdir", mode: "copy"

    input:
    path(database)
    val(outdir)
    //

    output:
    //

    script:
    """
    tandem.exe
    """
    //
}
