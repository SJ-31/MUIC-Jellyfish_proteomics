process ANNOTATE {
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_tsv)
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
