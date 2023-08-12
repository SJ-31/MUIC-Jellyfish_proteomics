process MSGF {
    publishDir "$outdir", mode: "copy"

    input:
    path(files)
    val(outdir)
    val(database)
    //
    output:
    path("${files.baseName}_msgf.tsv")
    //
    script: // -inst 3 specifies the Q-Exactive machine
    // You could change -tda to 1 to search the decoy database
    """
    cp $database ./database.fasta
    java -jar ~/tools/MSGFPlus/MSGFPlus.jar -s $files \
        -o temp.mzid \
        -d database.fasta \
        -inst 3 \
        -decoy rev \
        -t 20ppm \
        -minLength 7 \
        -m 3 \
        -addFeatures 1 \
        -maxMissedCleavages 2 \
        -tda 0
    Mzid_to_tsv_wrapper.sh temp.mzid ${files.baseName}_msgf.tsv
    rm database.fasta
    """
    //
}
/*
 * -tda 0 assumes that decoys are already present in the database */
/* -inst 1 specifies orbitrap */
/* -m 3 is HCD fragmentation
 */
