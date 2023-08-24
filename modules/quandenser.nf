process QUANDENSER {
    publishDir "$outdir", mode: "copy"
    errorStrategy 'ignore'

    input:
    path(mzMLs)
    //

    output:
    //

    shell:
    '''
    find $(pwd) -name "*.raw" \
            --filter "peakPicking true" \
        > peak_pick.txt
    quandenser -b peak_pick.txt -f .
    find . \\( -name "*.mzML" -o -name "*.mgf" \\) \
        -type f > temp_list.txt
    cat temp_list.txt | sed 's;./;!{outdir}/;' > !{params.pref}.falcon.txt

    '''
    //
}
