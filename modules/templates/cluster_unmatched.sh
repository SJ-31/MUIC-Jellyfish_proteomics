#!/usr/bin/env bash

for f in *.fasta; do
    cat "$f" >> tmp.txt
done
mv tmp.txt tmp.fasta
clean_fasta.py -i tmp.fasta -o unmatched.fasta

mkdir tmp
mmseqs createdb unmatched.fasta DB
mmseqs cluster DB clustered tmp --cov-mode 0
mmseqs createtsv DB DB clustered unmatched_clustered.tsv
mmseqs createseqfiledb DB clustered seq_clustered
mmseqs result2flat DB DB seq_clustered unmatched_clustered.fasta
