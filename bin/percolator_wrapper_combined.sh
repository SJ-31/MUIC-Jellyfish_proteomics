#!/usr/bin/env sh
# Wrapper for running percolator on a single combined pin file

while getopts "p:i:f:e:h" opt; do
    case $opt in
        p)
            prefix=$OPTARG ;;
        i)
            input=$OPTARG ;;
        f)
            fasta=$OPTARG ;;
        h)
            echo "  Usage"
            echo "  -p <prefix> -i <input_percolator_file>"
            echo "  -f <database_in_fasta_format>"
            echo "\n"
            ;;
        *) echo "Unsupported flag";;
    esac
done

percolator "$input" \
    -Y \
    -r "${prefix}"_percolator_peptides.tsv \
    -B "${prefix}"_percolator_decoy_peptides.tsv \
    -m "${prefix}"_percolator_psms.tsv \
    -M "${prefix}"_percolator_decoy_psms.tsv \
    -l "${prefix}"_percolator_proteins.tsv \
    -L "${prefix}"_percolator_decoy_proteins.tsv \
    --picked-protein "$fasta" \
    -P rev_ \
    --protein-enzyme trypsin \
    --protein-report-duplicates
cp .command.log "${prefix}"_percolator.log
#        -J "${prefix}"_features.tsv \
