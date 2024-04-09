process DEISOTOPE {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(mzMLs)
    val(outdir)
    //

    output:
    path("*mzML")
    path("*mgf")
    path("*msconvert.tsv")
    path("*log")
    //

    shell:
    def check = file("${outdir}/${params.pref}.msconvert.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/* .
        cp !{params.logdir}/deisotope.log .
        '''
    } else {
        '''
        mkdir files; mv *mzML files
        find $(pwd) -name "*.mzML" > convert.txt

        msconvert -f convert.txt \
                -o temp \
                --mzML \
                --filter "peakPicking true MS2Deisotope true MS2Denoise true"

        msconvert -f convert.txt \
                -o temp \
                --mgf \
                --filter "peakPicking true MS2Deisotope true MS2Denoise true"

        cd temp
        for file in *{mzML,mgf}; do
            mv $file msc_${file}
        done
        cd ..

        make_manifest.py -f temp -t "!{projectDir}/!{params.manifest_file}" \
            -o "!{params.pref}.msconvert.tsv" \
            -p "!{outdir}"

        mv temp/*{mzML,mgf} .
        cp .command.log deisotope.log
        '''
    }
}
