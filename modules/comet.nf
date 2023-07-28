process COMET {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzXMLs)
    val(outdir0
    //

    output:
    path ("${mzXML.baseName}_comet.tsv")
    //

    script:
    """
    cp $params.cometPars .
    philosopher workspace init
    comet --param default_comet.params $mzXMLs
    mv ${mzXMLs[0].baseName}.txt ${mzXMLs[0].baseName}_comet.tsv
    """
    // Comet does the reverse decoy search by default
}
