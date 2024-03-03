process COVERAGE_SPLIT {

    input:
    path(combined_tsv)
    //

    output:
    path("combined_slice*.tsv")
    //

    script:
    """
    Rscript ${params.bin}/R/protein_coverage.r \
        --split \
        --input $combined_tsv \
        --output_path .
    """
}

process COVERAGE_CALC {
    publishDir "${outdir}/.coverage", mode: "copy"

    input:
    path(sequences)
    val(outdir)
    //

    output:
    path("${sequences.baseName}_calculated.tsv")
    //

    script:
    def check = file("${outdir}/${sequences.baseName}_calculated.tsv")
    if (check.exists()) {
        """
        cp ${outdir}/${sequences.baseName}_calculated.tsv .
        """
    } else {
        """
        Rscript ${params.bin}/R/protein_coverage.r \
            --calculate \
            --input $sequences
        """
    }
    //
}

process COVERAGE_MERGE {
    publishDir "$outdir", mode: "copy"

    input:
    path(calculated_coverage)
    path(combined_file)
    val(outdir)
    //

    output:
    path("${combined_file.baseName}_wcoverage.tsv")
    //

    script:
    """
    Rscript ${params.bin}/R/protein_coverage.r \
        --merge \
        --input $combined_file \
        --input_path . \
        --output_path "${combined_file.baseName}_wcoverage".tsv \
        --alignment_file aligned_peptides.fasta
    """
}
