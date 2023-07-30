process MSFRAGGER {
    memory { 2.GB * task.attempt }
    publishDir "$outdir", mode: "copy"

    input:
    path(mzmls)
    val(pars)
    val(outdir)

    output:
    path("${params.pref}_msfragger.tsv")

    script:
    """
    java -Xmx32g -jar ~/tools/MSFragger-3.7/MSFragger-3.7.jar  \
        $pars \
        $mzmls
    mv ${params.pref}.tsv ${params.pref}_msfragger.tsv
    """
    // The output pin files need to be merged into a single one before doing the percolator analysis
    // Calibrate mass
}

// Something is wrong...
    // percolator_wrapper.sh $params.pref ${params.pref}.pin $params.databaseWdecoy msfragger
// for output: path("${params.pref}-msfragger_percolator")
