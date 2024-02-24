process COMBINE_ALL {
    publishDir "$outdir", mode: "copy"

    input:
    path(downloads)
    path(eggnog)
    path(interpro)
    path(directlfq)
    path(flashlfq)
    path(maxlfq)
    val(outdir)
    //

    output:
    path("${params.pref}_all.tsv"), emit: all
    path("aligned_peptides.fasta")
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
        --coverage TRUE \
        --sort_mods TRUE \
        --empai TRUE \
        --pfam2go $params.pfam2go \
        --interpro2go $params.interpro2go \
        --pfam_db $params.pfam_entries \
        --alignment_file aligned_peptides.fasta \
        --is_denovo $params.denovo \
        --fdr $params.fdr \
        --pep_thresh $params.pep_thresh \
        --output "${params.pref}_all.tsv" \
        --directlfq $directlfq \
        --maxlfq $maxlfq \
        --flashlfq $flashlfq \
        --r_source ${params.bin}/R \
        --python_source ${params.bin}
    """
    //
}
