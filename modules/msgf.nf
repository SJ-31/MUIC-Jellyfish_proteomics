process MSGF {
    publishDir "$outdir", mode: "copy"

    input:
    path(files)
    val(outdir)
    vald(database)
    //
    output:
    path("${files.baseName}_msgf.tsv")
    //
    script: // -inst 3 specifies the Q-Exactive machine
    // You could change -tda to 1 to search the decoy database
    """
    java -jar ~/tools/MSGFPlus/MSGFPlus.jar -s $files \
        -o temp.mzid \
        -d $database \
        -inst $params.inst \
        -decoy rev \
        -t ${params.masstolerance}ppm \
        -minLength $params.minpeplength \
        -m 3 \
        -addFeatures 1 \
        -maxMissedCleavages 2 \
        -tda 0
    Mzid_to_tsv_wrapper.sh temp.mzid ${files.baseName}_msgf.tsv
    """
    //
}
/*
 * -tda 0 assumes that decoys are already present in the database */
/* -inst 1 specifies orbitrap */
/* -m 3 is HCD fragmentation
 */
