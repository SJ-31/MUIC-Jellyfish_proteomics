#!/usr/bin/env sh

pin_header="SpecId	Label	ScanNr	Retention time	Charge	Score	Mass	Mass error [ppm]	m/z	Mass error [Da]	Isotope index	Combinatorics	Simple mass error [ppm]	PEP	Score	Precursor Intensity	Precursor apex fraction	Intensity coverage	Peak coverage	Mass deficit	Peptide	Proteins"

for file in *txt
do
    name=$(echo $file | sed 's/\..*//')
    format_mq.py -i $file -o $name.tsv
    sed 's/;/\t/g' $name.tsv > $name.pin
done

merge_tables.sh -r "${pin_header}" \
-o maxquant_all_pins.temp \
-p pin
