process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "{*.tsv,*.fasta}"
    publishDir "$outdir", mode: "copy", pattern: "${name}_eggnog-unmatched"
    publishDir "$outdir", mode: "copy", pattern: "${name}_interpro-unmatched"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("${params.pref}_downloads_anno-3.tsv"), emit: annotations
    path("${params.pref}_downloads_anno*")
    path("annotate_eggnog_unmatched")
    path("annotate_interpro_unmatched")
    path("*fasta"), optional: true
    //

    shell:
    template 'annotate.sh'
}
