process MS_MAPPING {
    publishDir "$outdir", mode: "copy"
    // Produces tsv file of MS/MS scans with precursor charge, retensionTime and
    //  intensity
    // Required by quantification scripts

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
    for file in *.mzML
        do
            Rscript !{params.bin}/ms_mapping.r $file ${file}_metrics.temp
    done

    merge_tables.sh -r "$metric_header" \
        -o  msms_scans.tsv \
        -p  temp
    '''
    //
}
