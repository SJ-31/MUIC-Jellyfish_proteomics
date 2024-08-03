process WRITE_QUANT {
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_results)
    path(aqreformat)
    val(outdir)
    //

    output:
    path("to_quantify.aq_reformat.tsv")
    //

    script:
    check = "to_quantify.aq_reformat.tsv"
    if (file("${outdir}/${check}").exists()) {
        """
        cp ${outdir}/to_quantify.aq_reformat.tsv .
        """
    } else {
    """
    helpers.py -t write_dlfq \
        -i $combined_results \
        -d $aqreformat \
        -o to_quantify.aq_reformat.tsv
    """
    }
    //
}

process LFQ_MERGE  {
    publishDir "$outdir", mode: "copy"

    input:
    tuple path(directlfq), path(top3), path(maxlfq)
    path(seq_header_mappings)
    val(outdir)
    //

    output:
    path("lfq_all.tsv")
    //

    script:
    check = "lfq_all.tsv"
    if (file("${outdir}/${check}").exists()) {
        """
        cp ${outdir}/lfq_all.tsv .
        """
    } else {
    """
    helpers.py -t merge \
        -d $directlfq \
        -p $top3 \
        -m $maxlfq \
        --seq_header_mapping $seq_header_mappings \
        -o lfq_all.tsv
    """
    }
    //
}
