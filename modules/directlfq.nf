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
    Rscript $params.bin/R/directlfq_format2.r \
        -p paths \
        -m $msms_mapping \
        -o directlfq.aq_reformat.tsv
    """
    //
}

process DIRECTLFQ {
    conda '/home/shannc/anaconda3/envs/directlfq'
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(aqreformat)
    val(outdir)
    val(logdir)
    //

    output:
    path("directlfq_prot.tsv"), emit: quant
    path("directlfq_ions.tsv")
    path("*.log")
    //

    script:
    def check = file("${outdir}/directlfq_prot.tsv")
    if (check.exists()) {
        """
        cp ${outdir}/directlfq_prot.tsv .
        cp ${outdir}/directlfq_ions.tsv .
        cp ${outdir}/directlfq.log .
        """
    } else {
        """
        directlfq lfq -i $aqreformat
        mv directlfq.aq_reformat.tsv.protein_intensities.tsv directlfq_prot.tsv
        mv directlfq.aq_reformat.tsv.ion_intensities.tsv directlfq_ions.tsv
        cp .command.out directlfq.log
        """
    }
    //
}
