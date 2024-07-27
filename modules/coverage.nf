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
    publishDir "$outdir", mode: "copy"

    input:
    path(sequences)
    path(seq_header_map)
    path(peptide_map)
    val(outdir)
    //

    output:
    tuple path(calculated), path(alignments)
    //

    script:
    calculated = "${sequences.baseName}_calculated.tsv"
    alignments = "${sequences.baseName}_alignments.tsv"
    if (file("${outdir}/${calculated}").exists()) {
        """
        cp "${outdir}/${alignments}" .
        cp "${outdir}/${calculated}" .
        """
    } else {
        """
        Rscript ${params.bin}/R/protein_coverage.r \
            --calculate \
            --input $sequences \
            --seq_header_map $seq_header_map \
            --peptide_map $peptide_map
        """
    }
}

process COVERAGE_MERGE {
    publishDir "$outdir", mode: "copy"

    input:
    path(calculated_coverage)
    path(combined_file)
    val(outdir)
    //

    output:
    path("${combined_file.baseName}_wcoverage.tsv"), emit: tsv
    path("aligned_peptides.fasta")
    path("aligned_peptides.tsv")
    path("alignment_metrics.tsv")
    path("all_mismatches.tsv")
    //

    script:
    """
    Rscript ${params.bin}/R/protein_coverage.r \
        --merge \
        --input $combined_file \
        --input_path . \
        --output_path "${combined_file.baseName}_wcoverage".tsv \
        --alignment_fasta aligned_peptides.fasta \
        --alignment_tsv aligned_peptides.tsv

    parse_alignment.py \
        -i aligned_peptides.fasta \
        -m alignment_metrics.tsv \
        -r all_mismatches.tsv
    """
}
