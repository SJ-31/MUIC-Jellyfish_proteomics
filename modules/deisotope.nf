process DEISOTOPE {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzMLs)
    val(outdir)
    //

    output:
    path("Proteowizard")
    path("*msconvert.txt")
    //

    shell:
    '''
    mkdir files; mv *mzML files
    find . -name "*.mzML" > convert.txt
    msconvert -f convert.txt \
            -o Proteowizard \
            --mzML \
            --filter "peakPicking true MS2Deisotope true MS2Denoise true"
    msconvert -f convert.txt \
            -o Proteowizard \
            --mgf \
            --filter "peakPicking true MS2Deisotope true MS2Denoise true"
    make_manifest.py -f Proteowizard -t !{projectDir}/!{params.manifest_file} \
    -o !{params.pref}.msconvert.txt \
        -p !{outdir}/Proteowizard
    '''
    //
}
