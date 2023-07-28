process MSGF {
    publishDir "$outdir", mode: "copy"

    input:
    path(files)
    val(outdir)
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
    -decoy $params.decoy_prefix
    -t $params.msgf_tolerance \
    -minLength $params.minpeplength \
    -m 3
    -addFeatures 1 \
    -maxMissedCleavages 2 \
    -tda 1
    Mzid_to_tsv_wrapper.sh temp.mzid ${files.baseName}_msgf.tsv
    """
    //
}
/*
 * -tda 1 means search a concatenated target-decoy database
 */
