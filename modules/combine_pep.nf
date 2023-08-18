process COMBINE_PEP {
    publishDir "$outdir", mode: "copy"
    tag "Combining PEPs: $percolator_out"

    input:
    path(percolator_out)
    val(is_psm)
    val(outdir)
    //

    output:
    path("$out")
    //
    script:
    psm_header = 'score q-value PEP peptide "Is_decoy"'
    protein_header = 'ProteinId ProteinGroupID q-value PEP Is_decoy'
    if (is_psm) {
        header = psm_header
        out = "combined_PEP_psms.tsv"
    } else {
        header = protein_header
        out = "combined_PEP_prot.tsv"
    }
    println "Combining PEPs: $percolator_out\n Output path: $outdir\n Output file: $out\n\n"
    """
    combine_pep_1_0_0.py \
        -i $percolator_out \
        -c $header \
        -o $out
    """
    //
}
