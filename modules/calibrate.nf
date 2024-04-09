process CALIBRATE {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"
    memory "10 GB"
    conda "/home/shannc/anaconda3/envs/metamorpheus"


    input:
    path(mzML) // Cannot be profile-mode mzML files
    path(database)
    val(outdir)
    //

    output:
    path("*mzML")
    path("*log")
    //

    shell:
    def check = file("${params.pref}.calibrated.txt")
    if (check.exists()) {
        '''
        cp !{outdir}/* .
        cp !{params.logdir}/calibrate.log .
        '''
    } else {
    '''
    mkdir to_calibrate
    mv *mzML to_calibrate
    mkdir output
    metamorpheus -s to_calibrate/* \
        -o output \
        -t !{params.config_dir}/CalibrateTaskconfig.toml \
        -d !{database}

    make_manifest.py -f output -t !{projectDir}/!{params.manifest_file} \
    -o !{params.pref}.calibrated.txt \
        -p !{outdir}

    mv output/*/*mzML .
    cp .command.log calibrate.log
    '''
    }
}
