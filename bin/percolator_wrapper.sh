#!/usr/bin/env sh

while getopts "p:i:f:e" opt; do
    case $opt in
        p)
            prefix=$OPTARG ;;
        i)
            input=$OPTARG ;;
        f)
            fasta=$OPTARG ;;
        e)
            engine=$OPTARG ;;
        *) echo "Unsupported flag";;
    esac
done
mkdir ${prefix}-${engine}_percolator
cd ${prefix}-${engine}_percolator
percolator ../"$input" \
    -r "${prefix}"_peptides.tsv \
    -B "${prefix}"_decoy_peptides.tsv \
    -m "${prefix}"_psms.tsv \
    -M "${prefix}"_decoy_psms.tsv \
    -l "${prefix}"_proteins.tsv \
    -L "${prefix}"_decoy_proteins.tsv \
    -J "${prefix}"_features.tsv \
    --picked-protein "$fasta" \
    --protein-decoy-pattern rev_ \
    --protein-enzyme trypsin \
    --protein-report-duplicates
