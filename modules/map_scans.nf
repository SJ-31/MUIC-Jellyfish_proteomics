process MAP_SCANS {
    publishDir "$outdir", mode: "copy"

    input:
    path(psm_file)
    val(mapping)
    val(outdir)
    //

    output:
    path("*scans.tsv")
    //

    shell:
    '''
    engine=$(echo !{psm_file.baseName} | sed 's/[-_].*//')
    Rscript !{params.bin}/get_scan_num.r \
        -i !{psm_file} \
        -m !{mapping}   \
        -e $engine \
        -o ${engine}_scans.tsv
    '''
    //
}
