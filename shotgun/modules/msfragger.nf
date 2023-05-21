process MSFRAGGER {
    publishDir "$outdir"

    input:
    path(files)
    val(outdir)

    output:
    path("${files.baseName}.pepXML")

    script:
    """
    java -jar ~/tools/MSFragger-3.7/MSFragger-3.7.jar  \
    ${projectDir}/config/MSFragger_params.params \
    $files
    """
}

