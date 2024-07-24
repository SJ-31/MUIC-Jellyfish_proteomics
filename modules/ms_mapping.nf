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
    output = "${outdir}/msms_scans.tsv"
    if (file(output).exists()) {
        '''
        cp !{output} .
        '''
    } else {
    '''
   header="scanNum	msLevel	retentionTime	precursorCharge	precursorIntensity	precursorMZ	totIonCurrent	peaksCount	basePeakMZ"

    for file in *.mzML
        do
            Rscript !{params.bin}/R/ms_mapping.r $file ${file}.temp
    done

    merge_tables.sh -r "$header" \
        -o  msms_scans.tsv \
        -p  temp
    '''
    }
    //
}
