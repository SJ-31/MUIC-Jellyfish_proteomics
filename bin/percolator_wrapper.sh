#!/usr/bin/env sh

while getopts "p:i:f:e:h" opt; do
    case $opt in
        p)
            prefix=$OPTARG ;;
        i)
            input=$OPTARG ;;
        f)
            fasta=$OPTARG ;;
        e)
            engine=$OPTARG ;;
        h)
            echo "  Usage"
            echo "  -p <prefix> -i <input_percolator_file>"
            echo "  -f <database_in_fasta_format> -e <engine>"
            echo "\n"
            ;;
        *) echo "Unsupported flag";;
    esac
done

mkdir ${prefix}-${engine}_percolator

merge_tables () {
    echo -e "$1" > "$2"
    find . -maxdepth 1 -regex "[^_]*_${3}.tsv" \
        ! -name "${engine}"_"${3}".tsv \
        -exec sed -se 1d {} + >> \
        "$2"
}

for i in *$input
do
    name=$(echo $i | sed 's/\..*//')
    percolator "$i" \
        -Y \
        -r "${name}"_peptides.tsv \
        -B "${name}"_decoy_peptides.tsv \
        -m "${name}"_psms.tsv \
        -M "${name}"_decoy_psms.tsv \
        -l "${name}"_proteins.tsv \
        -L "${name}"_decoy_proteins.tsv \
        --picked-protein "$fasta" \
        -P rev_ \
        --protein-enzyme trypsin \
        # --protein-report-duplicates
#        -J "${name}"_features.tsv \
done

find . -maxdepth 1 -regex '\./.*\(peptides\|proteins\|features\|psms\).*' -exec mv {} ${prefix}-${engine}_percolator/ \;

peptide_header="PSMId	score	q-value	posterior_error_prob	peptide	proteinIds"
feature_header="SpecId	Label	ScanNr	ExpMass	CalcMass	rank	abs_ppm	isotope_errors	log10_evalue	hyperscore	delta_hyperscore	matched_ion_num	complementary_ions	ion_series	weighted_average_abs_fragment_ppm	peptide_length	ntt	nmc	charge_1	charge_2	charge_3	charge_4	charge_5	charge_6	charge_7_or_more	15.994915M	Peptide	Proteins"
protein_header="ProteinId	ProteinGroupId	q-value	posterior_error_prob	peptideIds"

cd ${prefix}-${engine}_percolator
merge_tables "$peptide_header" "${engine}"_peptides.tsv peptides
merge_tables "$peptide_header" "${engine}"_decoy_peptides.tsv decoy_peptides
merge_tables "$peptide_header" "${engine}"_psms.tsv psms
merge_tables "$peptide_header" "${engine}"_decoy_psms.tsv decoy_psms
merge_tables "$protein_header" "${engine}"_proteins.tsv proteins
merge_tables "$protein_header" "${engine}"_decoy_proteins.tsv decoy_proteins
mv *${engine}* ..
# merge_tables "$feature_header" "${prefix}"_features.tsv features
