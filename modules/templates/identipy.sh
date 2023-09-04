pin_header="ScanID	Label	ScanNr	fragmentMT	sumI	Expect	Hyperscore	X..proteins	Missed.cleavages	Mass.difference	Calculated.mass	Total.ions	Matched.ions	Rank	compensation_voltage	RT	Assumed.charge	Peptide	Proteins"
tsv_header="Title	Assumed charge	RT	compensation_voltage	Rank	Matched ions	Total ions	Calculated mass	Mass difference	Missed cleavages	Proteins	# proteins	Sequence	Modified sequence	Hyperscore	Expect	sumI	fragmentMT"
database="!{database}"
cp !{params.config}/identipy.cfg .
cat identipy.cfg | sed "s;^database:.*;database: $database;" > config.cfg

for mzML in *mzML
do
    identipy $mzML \
        -cfg config.cfg \
        -at \
        -out .
done

for pep in *.pep.xml
do
    name=$(echo $pep | sed 's/\..*//')
    identipy2pin $pep
    mv ${name}.pin ${name}.temp
    awk -v n="$name" '{
        if (FNR == 1)
            {print $0; next}
        else
            {$1 = name"."$1; print $0}
        }' ${name}.temp > ${name}.pin
done

merge_tables.sh -r "$pin_header" \
    -o identipy_all_pins.temp \
    -p pin
