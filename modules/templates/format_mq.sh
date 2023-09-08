#!/usr/bin/env sh
# pin_header="SpecId	Label	ScanNr	Retention time	Charge	Score	Mass	m/z	Score diff	Precursor Intensity	Length	Intensity coverage	Peak coverage	Number of matches	Delta score	Missed cleavages	Precursor apex fraction	Mass error [ppm]	Peptide	Proteins"
format_mq.py -i !{maxquant_msms} -o temp.pin
sed 's/;/\t/g' temp.pin > maxquant_all_pins.temp

# for file in *txt
# do
#     name=$(echo $file | sed 's/\..*//')
#     format_mq.py -i $file -o $name.pin
#     sed 's/;/\t/g' $name.tsv > $name.pin
# done

# merge_tables.sh -r "${pin_header}" \
#     -o maxquant_all_pins.temp \
#     -p pin
