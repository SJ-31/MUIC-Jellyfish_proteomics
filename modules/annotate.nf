process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "*.tsv"
    errorStrategy 'ignore'

    input:
    path(combined_tsv)
    path(map_file)
    val(outdir)
    //

    output:
    path("annotated_proteins.tsv"), emit: annotations
    path("denovo_matches.fasta"), emit: denovo
    path("transcriptome_matches.fasta"), emit: transcriptome
    //

    script:
    """
    annotate.py $combined_tsv annotated_proteins.tsv $map_file
    """
    //
}
