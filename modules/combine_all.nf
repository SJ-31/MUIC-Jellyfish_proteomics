process COMBINE_ALL {
    publishDir "$outdir", mode: "copy", pattern: "*.tsv"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(downloads)
    path(eggnog)
    path(interpro)
    path(directlfq)
    path(flashlfq)
    path(maxlfq)
    val(outdir)
    val(logdir)
    //

    output:
    path("${params.pref}_all.tsv"), emit: all
    path("combine_all.log")
    //

    script:
    if (eggnog == "no_file") {
        eggnog = interpro = ""
    }
    """
    Rscript ${params.bin}/R/combine_all.r \
        --eggnog $eggnog \
        --interpro $interpro \
        --downloads $downloads \
        --sort_mods TRUE \
        --empai FALSE \
        --pfam2go $params.pfam2go \
        --interpro2go $params.interpro2go \
        --kegg2go $params.kegg2go \
        --ec2go $params.ec2go \
        --pfam_db $params.pfam_entries \
        --is_denovo $params.denovo \
        --fdr $params.fdr \
        --pep_thresh $params.pep_thresh \
        --output "${params.pref}_all.tsv" \
        --directlfq $directlfq \
        --maxlfq $maxlfq \
        --flashlfq $flashlfq \
        --r_source ${params.bin}/R \
        --python_source ${params.bin}

    cat .command.log >> combine_all.log
    """
    //
}
