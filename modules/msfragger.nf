process MSFRAGGER {
    memory { 2.GB * task.attempt }
    publishDir "$outdir", mode: "copy"

    input:
    path(files)
    val(outdir)

    output:
    path("${files.baseName}_msfragger.tsv")

    script:
    """
    java -jar ~/tools/MSFragger-3.7/MSFragger-3.7.jar  \
    ${projectDir}/config/MSFragger_params.params \
    $files
    mv ${files.baseName}.tsv ${files.baseName}_msfragger.tsv
    """
}

