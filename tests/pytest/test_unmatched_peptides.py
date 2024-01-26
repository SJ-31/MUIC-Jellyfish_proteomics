#!/usr/bin/env ipython
import os
import sys

workdir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
output = f"{workdir}/tests/pytest/output"
m = "~/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish"
import unmatched_peptides as up

args = {
    "pep_threshold": 1,
    "q_threshold": 0.05,
    "output": "unmatched_peptides.fasta",
    "input_path": f"{workdir}/results/jellyfish/1-First_pass/Percolator",
    "unmatched_tsv": "./results/jellyfish/1-First_pass/Unmatched/unmatched_peptides.tsv",
}

up.main(args)
result = pd.read_csv(args["unmatched_tsv"], sep="\t")
