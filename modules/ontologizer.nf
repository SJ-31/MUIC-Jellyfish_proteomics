process ONTOLOGIZER {
    // For overrepresentation analysis
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("*")
    path("ontologizer-*"), emit: over
    //

    script:
    """
    Rscript $params.bin/R/ontologizer.r \
        -i $combined_tsv

    java -jar $params.ontologizer_jar -g $params.go_ontology \
            -a protein_mappings.ids \
            -p universe.txt \
            -s id_with_open.txt \
            -m "Bonferroni-Holm"

    java -jar $params.ontologizer_jar -g $params.go_ontology \
            -a protein_mappings.ids \
            -p universe.txt \
            -s unknown_to_db.txt \
            -m "Bonferroni-Holm"
    mv table-unknown_to_db* ontologizer-unknown_to_db.txt
    mv table-id_with_open* ontologizer-id_with_open.txt
    """
    //
}
