process COMET {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(pref), path(mzXMLs)
    val(outdir)
    //

    output:
    path ("${params.pref}_comet.tsv")
    path("${params.pref}-comet_percolator")
    //

    script:
    """
    cp $params.cometPars .
    philosopher workspace init
    comet --param default_comet.params $mzXMLs
    percolator_wrapper.sh $params.pref ${params.pref}.pin $params.databaseWdecoy comet
    mv ${params.pref}.txt ${params.pref}_comet.tsv
    """
    // Comet does the reverse decoy search by default
}
