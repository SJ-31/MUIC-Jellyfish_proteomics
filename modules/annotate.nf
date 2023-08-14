process ANNOTATE {
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("annotated_proteins.tsv")
    //

    script:
    """
    annotate.py $combined_tsv annotated_proteins.tsv
    """
    //
}
