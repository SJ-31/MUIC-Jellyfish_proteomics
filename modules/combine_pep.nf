process COMBINE_PEP {
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_PEPs)
    val(is_psm)
    //

    output:
    //

    script:
    peptide_header = 'score q-value PEP peptide "Is decoy"'
    protein_header = 'ProteinId ProteinGroupdId q-value PEP "Is decoy"'
    if (is_peptide) {
        header = peptide_header
    }
    """
    combine_pep_1_0_0.py \
        -i \
        -c $header \
        -e $engines
    """
    //
}
