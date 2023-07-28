#!/bin/bash
cat $@ | seqkit rmdup  > inter1.fasta
if [[ -e combined_decoys.fasta ]]
then
    rm combined_decoys.fasta
fi
$BIN/create_decoys.py inter1.fasta combined_decoys.fasta header_mappings.tsv
rm inter1.fasta
