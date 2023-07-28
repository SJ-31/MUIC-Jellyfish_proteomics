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
        $pars \
        $files
    mv ${files.baseName}.tsv ${files.baseName}_msfragger.tsv
    """
    // Calibrate mass
}

