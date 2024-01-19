tsv_header="Title	Assumed charge	RT	compensation_voltage	Rank	Matched ions	Total ions	Calculated mass	Mass difference	Missed cleavages	Proteins	# proteins	Sequence	Modified sequence	Hyperscore	Expect	sumI	fragmentMT"
database="!{database}"
cp !{params.config_dir}/identipy.cfg .
cat identipy.cfg | sed "s;^database:.*;database: $database;" > config.cfg

for mzML in *mzML
do
    name=$(echo "$mzML" | sed 's/.mzML/.pep.xml/' )
    if [ -e !{outdir}/"$name" ]; then
        echo "${mzML} already completed!"
        cp !{outdir}/${name} .
    else
    identipy ${mzML} \
        -cfg config.cfg \
        -at \
        -out .
    fi
done
