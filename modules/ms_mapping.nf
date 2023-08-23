process MS_MAPPING {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzmls)
    val(outdir)
    //

    output:
    path("msms_scans.tsv")
    //

    shell:
    '''
    metric_header="scanNum	retensionTime	precursorCharge	precursorIntensity"
    for file in *.mzml
        do
            Rscript !{params.bin}/ms_metrics.r $file ${file}_metrics.temp
    done

    merge_tables.sh -r "$metric_header" \
        -o  msms_scans.tsv \
        -p  temp
    '''
    //
}
