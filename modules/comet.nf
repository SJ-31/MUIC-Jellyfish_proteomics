process COMET {
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(indexed_mzMLs)
    val(outdir)
    val(logdir)
    val(database)
    //

    output:
    path ("${params.pref}_comet.tsv")
    tuple val("comet"), path("comet_all_pins.temp"), emit: percolator
    path("*.log")
    //

    shell:
    def check = file("${outdir}/${params.pref}_comet.tsv")
    if (check.exists()) {
        '''
        cp !{outdir}/!{params.pref}_comet.tsv .
        cp !{outdir}/comet_all_pins.temp .
        cp !{logdir}/comet.log .
        '''
    } else {
        '''
        tsv_header="scan	num	charge	exp_neutral_mass	calc_neutral_mass	e-value	xcorr	delta_cn	sp_score	ions_matched	ions_total	plain_peptide	modified_peptide	prev_aa	next_aa	protein	protein_count	modifications"
        pin_header="SpecId	Label	ScanNr	ExpMass	CalcMass	lnrSp	deltLCn	deltCn	lnExpect	Xcorr	Sp	IonFrac	Mass	PepLen	Charge1	Charge2	Charge3	Charge4	Charge5	Charge6	enzN	enzC	enzInt	lnNumSP	dM	absdM	Peptide	Proteins"
        database="!{database}"

        cp !{params.config}/default_comet.params .
        cat default_comet.params | sed "s;database.*;database_name = $database;" > comet.params
        philosopher workspace --init
        philosopher comet --param comet.params !{indexed_mzMLs} > comet.log
        echo "$tsv_header" > "!{params.pref}_comet.tsv"
        find . -maxdepth 1 -name "*txt" -exec sed -se 2d {} + >> \
            "!{params.pref}_comet.tsv"

        merge_tables.sh -r "$pin_header" \
            -o comet_all_pins.temp \
            -p pin
        '''
    }
    //
}
