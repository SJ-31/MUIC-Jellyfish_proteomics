#!/usr/bin/env bash

path="../data/protein_databases/denovo_normal_mgf"
cat "${path}/*.fasta" > all_denovo.fasta

../bin/helpers.py -t denovo_stats -i all_denovo.fasta -s ../results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv
