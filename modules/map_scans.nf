process MAP_SCANS {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(engine), path(psm_file), path(proteins)
    path(unmatched_peptides)
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
        -p !{proteins} \
        -u !{unmatched_peptides} \
        -m !{mapping}   \
        -e $engine \
        -o ${engine}_scans.tsv
    '''
    //
}
