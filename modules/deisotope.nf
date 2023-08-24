process DEISOTOPE {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzMLs)
    path(outdir)
    //

    output:
    path("*")
    //

    shell:
    '''
    mkdir mzML; mv *mzML mzML
    find $(pwd)/mzML -name "*.mzML" > convert.txt
    for format in {mzML,mgf}
    do
        msconvert -f convert.txt \
            -o . \
            --"${format}" \
            --filter "peakPicking true MS2Deisotope true MS2Denoise true"
    done
    find . \\( -name "*.mzML" -o -name "*.mgf" \\) \
        -type f > temp_list.txt
    cat temp_list.txt | sed 's;./;!{outdir}/;' > !{params.pref}.msconvert.txt
    '''
    //
}
