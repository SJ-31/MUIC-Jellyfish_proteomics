process CALIBRATE {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logs", mode: "copy", pattern: "${params.pref}_results*"
    memory "10 GB"
    conda "/home/shannc/anaconda3/envs/metamorpheus"


    input:
    path(raw)
    path(database)
    val(outdir)
    //

    output:
    path("*mzML")
    //

    script:
    """
    metamorpheus -s *raw \
        -o . \
        -t ${params.config}/CalibrateTaskconfig.toml \
        -d $database

    find . -name "*mzML" -exec mv {} . \;
    """
    //
}
