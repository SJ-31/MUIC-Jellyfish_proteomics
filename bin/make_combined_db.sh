#!/bin/bash
cat $@ | seqkit rmdup > inter1.fasta
if [[ -e decoysWnormal.fasta ]]
then
    rm decoysWnormal.fasta all_decoys.fasta all_normal.fasta
fi
$BIN/create_decoys.py inter1.fasta seq-header_mappings.tsv
rm inter1.fasta
