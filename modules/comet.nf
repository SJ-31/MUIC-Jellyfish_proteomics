process COMET {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzXMLs)
    val(outdir)
    //

    output:
    path ("${params.pref}_comet.tsv")
    path ("comet*.tsv")
    //

    script:
    tsv_header = "scan	num	charge	exp_neutral_mass	calc_neutral_mass	e-value	xcorr	delta_cn	sp_score	ions_matched	ions_total	plain_peptide	modified_peptide	prev_aa	next_aa	protein	protein_count	modifications"
    pin_header = "SpecId	Label	ScanNr	ExpMass	CalcMass	lnrSp	deltLCn	deltCn	lnExpect	Xcorr	Sp	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	Charge4	Charge5	Charge6	enzN	enzC	enzInt	lnNumSP	dM	absdM	Peptide	Proteins"
    """
    cp $params.cometPars .
    philosopher workspace --init
    philosopher comet --param default_comet.params $mzXMLs

    echo -e "$tsv_header" > "${params.pref}_comet.tsv"
    find . -maxdepth 1 -name "*txt" -exec sed -se 2d {} + >> \
        "${params.pref}_comet.tsv"

    merge_tables.sh -r "$pin_header" \
        -o comet_all_pins.temp \
        -p pin

    percolator_wrapper_combined.sh \
        -p comet \
        -i comet_all_pins.temp \
        -f $params.database \
        -e comet

    Rscript $params.bin/to_combined_PEP.r \
        -m comet_percolator_psms.tsv \
        -d comet_percolator_decoy_psms.tsv \
        -o comet_psm2combined_PEP.tsv

    Rscript $params.bin/to_combined_PEP.r \
        -m comet_percolator_proteins.tsv \
        -d comet_percolator_decoy_proteins.tsv \
        --protein_matches \
        -o comet_prot2combined_PEP.tsv
    """
    //
}
