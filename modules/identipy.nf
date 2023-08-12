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
    tsv_header="Title	Assumed charge	RT	compensation_voltage	Rank	Matched ions	Total ions	Calculated mass	Mass difference	Missed cleavages	Proteins	# proteins	Sequence	Modified sequence	Hyperscore	Expect	sumI	fragmentMT"
    database="!{database}"
    cp !{params.config}/identipy.cfg .
    cat identipy.cfg | sed "s;^database:.*;database: $database;" > config.cfg
    identipy_wrapper.sh -g config.cfg \
        -p !{params.pref} \
        -c !{params.bin}/identipy2pin.r

    merge_tables.sh -r "$pin_header" \
        -o identipy_all_pins.temp \
        -p pin

    merge_tables.sh -r "$tsv_header" \
        -o "!{params.pref}"_identipy.txt \
        -p tsv
    mv "!{params.pref}"_identipy.txt "!{params.pref}"_identipy.tsv
    '''
    // -at Auto-tune search parameters
    // -mc Number of missed cleavages
    // -method {reverse,shuffle} for decoy generation
    // -prefix decoy prefix
}
