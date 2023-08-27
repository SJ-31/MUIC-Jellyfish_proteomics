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
    for format in {mzML,mgf}
    do
        msconvert -f convert.txt \
            -o Proteowizard \
            --"${format}" \
            --filter "peakPicking true MS2Deisotope true MS2Denoise true"
    done
    make_manifest.py -f . -t !{params.manifest_file} \
    -o !{params.pref}.msconvert.txt
    '''
    //
}
