process MSGF {
    input:
    path(files)
    //
    output:
    
    //
    script: // -inst 3 specifies the Q-Exactive machine
    // You could change -tda to 1 to search the decoy database
    """
    java -jar ~/tools/MSGFPlus/MSGFPlus.jar -s $files \
    -o <put output file here> \
    -d $params.database \
    -t <precursor mass tolerance> \
    -inst 3 -tda 0

    """
    //
}
