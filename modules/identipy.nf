process IDENTIPY {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/identipy"

    input:
    path(mzMLs)
    val(outdir)
    val(database)
    //
    output:
    path("${params.pref}_identipy.tsv")
    tuple val("identipy"), path("identipy_all_pins.temp"), emit: percolator
    //
    shell:
    '''
    pin_header="ScanID	Label	ScanNr	fragmentMT	sumI	Expect	Hyperscore	X..proteins	Missed.cleavages	Mass.difference	Calculated.mass	Total.ions	Matched.ions	Rank	compensation_voltage	RT	Assumed.charge	Peptide	Proteins"
    database="!{database}"
    cp !{params.config/identipy.cfg} .
    cat identipy.cfg | sed "s/^database:.*/database: $database/" > config.cfg
    identipy_wrapper.sh -g config.cfg \
        -p $params.pref \
        -c $params.bin/identipy2pin.r

    merge_tables.sh -r "$pin_header" \
        -o identipy_all_pins.temp \
        -p pin
    '''
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -method {reverse,shuffle} for decoy generation
    // -prefix decoy prefix
}
