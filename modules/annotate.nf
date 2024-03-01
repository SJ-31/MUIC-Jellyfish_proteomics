process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "{*.tsv,*.fasta}"
    publishDir "$outdir", mode: "copy", pattern: "annotate_eggnog-unmatched"
    publishDir "$outdir", mode: "copy", pattern: "annotate_interpro-unmatched"

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
    check = file("${outdir}/${params.pref}_downloads_anno-3.tsv")
    if (check.exists()) {
        '''
        cp -r !{outdir}/* .
        '''
    } else {
        template 'annotate.sh'
    }
}
