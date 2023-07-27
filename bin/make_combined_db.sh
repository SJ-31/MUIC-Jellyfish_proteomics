#!/bin/bash
cat $@ | seqkit rmdup  > inter1.fasta
$BIN/header_mapping.sh inter1.fasta > decoy_headings.tsv # Need this intermediate because the 
# Uniprot intermediates are removed in the next step
cat inter1.fasta | sed 's/>.*|\(.*\)|.*/>\1/' > inter2.fasta
$BIN/simplify_fasta.sh inter2.fasta > combined.fasta
if [[ -e combined_decoys.fasta ]]
then
    rm combined_decoys.fasta
fi
$BIN/create_decoys.py combined.fasta decoys.fasta
# $BIN/make_decoy_ursgal.py combined.fasta  combined_decoys.fasta
cat combined.fasta decoys.fasta > combined_decoys.fasta
rm inter1.fasta inter2.fasta decoys.fasta
