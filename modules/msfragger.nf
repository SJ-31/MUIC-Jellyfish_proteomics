process MSFRAGGER {
    memory { 2.GB * task.attempt }
    publishDir "$outdir", mode: "copy"

    input:
    path(files)
    val(pars)
    val(outdir)

    output:
    path("${files.baseName}_msfragger.tsv")

    script:
    """
    java -Xmx32g -jar ~/tools/MSFragger-3.7/MSFragger-3.7.jar  \
        --database_name $params.databaseWdecoy \
        $pars \
        ${projectDir}/config/MSFragger_params.params \
        $files
    mv ${files.baseName}.tsv ${files.baseName}_msfragger.tsv
    """
    // Calibrate mass
}

