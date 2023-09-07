process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "*.tsv"
    errorStrategy 'ignore'

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("annotated_proteins.tsv"), emit: annotations
    //

    script:
    """
    annotate.py -i $combined_tsv \
    -o annotated_proteins.tsv \
    -s $map_file
    """
    //
}
