#!/usr/bin/env sh
while getopts "d:p:c:" opt; do
    case $opt in
        d)
            database=$OPTARG ;;
        p)
            prefix=$OPTARG ;;
        c)
            converter=$OPTARG ;;
    esac
done
#
tsv_header="Title	Assumed charge	RT	compensation_voltage	Rank	Matched ions	Total ions	Calculated mass	Mass difference	Missed cleavages	Proteins	# proteins	Sequence	Modified sequence	Hyperscore	Expect	sumI	fragmentMT"
for mzML in *mzML
do
    identipy $mzmL \
        -db "$database"
        -of tsv \
        -at \
        -prefix rev_
        -out .
done

merge_tables.sh -r "$tsv_header" \
    -o "${prefix}"_identipy.txt \
    -p tsv
mv "${prefix}"_identipy.txt "${prefix}"_identipy.tsv

for tsv in *tsv
do
    sample=$(echo $tsv | sed 's/\.tsv//')
     Rscript "$converter" $tsv \
    ${sample}.pin \
    ${sample}
done
