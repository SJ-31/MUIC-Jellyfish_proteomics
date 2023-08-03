process COMET {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzXMLs)
    val(outdir)
    //

    output:
    path ("${params.pref}_comet.tsv")
    path ("comet*.tsv")
    tuple val("comet"), path("comet_all_pins.temp"), emit: percolator
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
    """
    //
}
