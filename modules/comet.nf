process COMET {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzXML)
    val(outdir0
    //

    output:
    path ("${mzXML.baseName}_comet.tsv")
    //

    script:
    """
    cp $params.cometPars .
    philosopher workspace init
    comet --param default_comet.params $mzXML
    mv ${mzXML.baseName}.txt ${mzXML.baseName}_comet.tsv
    """
    // Comet does the reverse decoy search by default
}
