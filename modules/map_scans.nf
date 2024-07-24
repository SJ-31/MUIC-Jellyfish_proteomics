process MAP_SCANS {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(engine), path(psm_file), path(proteins)
    val(mapping)
    val(outdir)
    //

    output:
    path("*scans.tsv")
    //

    shell:
    check = "${outdir}/${engine}_scans.tsv"
    if (file(check).exists()) {
        '''
        cp "!{check}" .
        '''
    } else {
        '''
        engine=$(echo !{psm_file.baseName} | sed 's/[-_].*//')
        Rscript !{params.bin}/R/get_scan_num.r \
            -i !{psm_file} \
            -p !{proteins} \
            -m !{mapping}   \
            -e $engine \
            -o ${engine}_scans.tsv
        '''
    }
    //
}
