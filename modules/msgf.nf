process MSGF {
    publishDir "$outdir", mode: "copy"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(file)
    val(outdir)
    val(logdir)
    val(database)
    //
    output:
    path("${params.pref}_msgf.pin").out.pin
    path("*.log")
    //
    shell: // -inst 3 specifies the Q-Exactive machine
    // You could change -tda to 1 to search the decoy database
    def check = file("${outdir}/${params.pref}_msgf.pin")
    if (check.exists()) {
        '''
        cp !{outdir}/* .
        cp !{logdir}/msgf* .
        '''
    } else {
        '''
        cp !{database} wdecoys.fasta
        seqkit grep wdecoys.fasta -n -r -p "rev_" -v > normal.fasta
        seqkit grep wdecoys.fasta -n -r -p "rev_" > decoys.fasta
        names=("decoys" "normal")
        for n in ${names[@]}; do
            java -jar !{params.msgf} -s $files \
                -o ${n}.mzid \
                -d ${n}.fasta \
                -inst 3 \
                -decoy rev \
                -t 20ppm \
                -minLength 7 \
                -m 3 \
                -addFeatures 1 \
                -maxMissedCleavages 2 \
            -tda 0 > msgf_${n}_search.log
        done
        rm *.fasta
        msgf2pin.py \
            -d decoys.mzid \
            -v normal.mzid \
            -o !{params.pref}_msgf.pin
        '''
    }
}
/*
 * -tda 0 means don't create new decoys
/* -inst 1 specifies orbitrap */
/* -m 3 is HCD fragmentation
 */

process MSGF_MERGE {
    publishDir "$outdir", mode: "copy"

    input:
    path(pin_files)
    val(outdir)
    //

    output:
    tuple val("msgf"), path("msgf_all_pins.temp")
    //

    shell:
    '''
    pin_header="SpecId	Label	ScanNr	ExpMass	CalcMass	RawScore	DeNovoScore	ScoreRatio	Energy	lnEValue	IsotopeError	lnExplainedIonCurrentRatio	lnNTermIonCurrentRatio	lnCTermIonCurrentRatio	lnMS2IonCurrent	Mass	PepLen	dM	absdM	MeanErrorTop7	sqMeanErrorTop7	StdevErrorTop7	Charge1	Charge2	Charge3	Charge4	Charge5	enzN	enzC	enzInt	Peptide	Proteins"
    merge_tables.sh -r "$pin_header" \
        -o msgf_all_pins.temp \
        -p pin
    '''
    //
}
