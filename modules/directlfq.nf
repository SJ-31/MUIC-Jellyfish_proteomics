process DIRECTLFQ_FORMAT {
    publishDir "$outdir", mode: "copy"

    input:
    path(scan_prot_mappings)
    path(msms_mapping)
    val(outdir)
    //

    output:
    path("directlfq.aq_reformat.tsv")
    //

    script:
    """
    mkdir paths; mv $scan_prot_mappings paths
    Rscript $params.bin/directlfq_format2.r \
        -p paths
        -o directlfq.aq_reformat.tsv
    """
    //
}

process DIRECTLFQ {
    conda '/home/shannc/anaconda3/envs/directlfq'
    publishDir "$outdir", mode: "copy"

    input:
    path(aqreformat)
    val(outdir)
    //

    output:
    path("directlfq_prot.tsv")
    //

    script:
    """
    directlfq lfq -i $aqreformat
    mv directlfq.aq_reformat.tsv.protein_intensities.tsv directlfq_prot.tsv
    """
    //
}
