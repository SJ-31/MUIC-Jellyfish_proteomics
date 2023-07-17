process MSGF {
    input:
    path(files)
    //
    output:
    path("${files.baseName}_msgf.tsv")
    //
    script: // -inst 3 specifies the Q-Exactive machine
    // You could change -tda to 1 to search the decoy database
    """
    java -jar ~/tools/MSGFPlus/MSGFPlus.jar -s $files \
        -o temp.mzid \
    -d $params.database \
    -inst $params.inst \
    -t $params.msgf_tolerance \
    -tda 0
    mono Mzid_to_tsv_wrapper.sh temp.mzid ${files.baseName_msgf.tsv}
    """
    //
}
