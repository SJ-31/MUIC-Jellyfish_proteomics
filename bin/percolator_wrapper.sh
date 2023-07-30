#!/usr/bin/env sh

prefix=$1
input=$2
fasta=$3
engine=$4
percolator "$input" \
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
mkdir ${prefix}_percolator
find . -maxdepth 1 -type f \! -name "$input" -exec mv {} ${prefix}-${engine}_percolator \;
