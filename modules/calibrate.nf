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
    path("*mgf")
    path("*log")
    path("${params.pref}.calibrated.tsv")
    //

    shell:
    def check = file("${outdir}/${params.pref}.calibrated.tsv")
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

    mv output/*/*mzML .

    find $(pwd)/ -name "*-calib.mzML" > convert.txt

    msconvert -f convert.txt \
            -o . \
            --mgf

    make_manifest.py -f . -t !{projectDir}/!{params.manifest_file} \
    -o !{params.pref}.calibrated.tsv \
        -p !{outdir}

    cp .command.log calibrate.log
    '''
    }
}
