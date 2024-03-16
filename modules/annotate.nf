process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "{*.tsv,*.fasta}"
    publishDir "$outdir", mode: "copy", pattern: "annotate_eggnog_unmatched"
    publishDir "$outdir", mode: "copy", pattern: "annotate_interpro_unmatched"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("${params.pref}_downloads_anno-3.tsv"), emit: annotations
    path("${params.pref}_downloads_anno*")
    path("annotate_eggnog_unmatched"), optional: true
    path("annotation_complete.fasta"), emit: complete, optional: true
    path("annotate_interpro_unmatched"), optional: true
    path("still_unannotated.fasta"), emit: unannotated, optional: true
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
