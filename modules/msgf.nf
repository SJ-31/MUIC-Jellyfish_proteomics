process MSGF {
    publishDir "$outdir", mode: "copy", pattern: "{*.pin,*.mzid}"
    publishDir "$logdir", mode: "copy", pattern: "*.log"

    input:
    path(file)
    val(outdir)
    val(logdir)
    val(database)
    //
    output:
    path("${params.pref}-${file.baseName}_msgf.pin"), emit: pin
    path("*.mzid")
    path("*.log")
    //
    shell: // -inst 3 specifies the Q-Exactive machine
    // You could change -tda to 1 to search the decoy database
    name = "${params.pref}-${file.baseName}"
    check = file("${outdir}/${name}_msgf.pin")
    if (check.exists()) {
        '''
        cp !{outdir}/*!{file.baseName}* .
        cp !{logdir}/*!{file.baseName}* .
        '''
    } else {
        '''
        cp !{database} wdecoys.fasta
        seqkit grep wdecoys.fasta -n -r -p "rev_" -v > normal.fasta
        seqkit grep wdecoys.fasta -n -r -p "rev_" > decoys.fasta
        names=("decoys" "normal")
        for n in ${names[@]}; do
            java -jar !{params.msgf} -s !{file} \
                -o ${n}-!{file.baseName}.mzid \
                -d ${n}.fasta \
                -conf !{params.config_dir}/MSGFPlus_Params.txt \
                -addFeatures 1 > !{name}_search.log
        done
        rm *.fasta
        msgf2pin.py \
            -d decoys-!{file.baseName}.mzid \
            -v normal-!{file.baseName}.mzid \
            -o !{name}_msgf.pin
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
    output=msgf_all_pins.temp
    pin_header="SpecId	Label	ScanNr	ExpMass	CalcMass	RawScore	DeNovoScore	ScoreRatio	Energy	lnEValue	IsotopeError	lnExplainedIonCurrentRatio	lnNTermIonCurrentRatio	lnCTermIonCurrentRatio	lnMS2IonCurrent	Mass	PepLen	dM	absdM	MeanErrorTop7	sqMeanErrorTop7	StdevErrorTop7	Charge1	Charge2	Charge3	Charge4	Charge5	enzN	enzC	enzInt	Peptide	Proteins"
    echo "$pin_header" > "$output"
    find . -maxdepth 1 -name "*pin" -exec sed -se 1,2d {} + >> "$output"
    '''
    //
}
