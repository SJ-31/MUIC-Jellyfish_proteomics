#!/usr/bin/env python
import pandas as pd
import sys

sys.path.append("./bin")
import sort_blast as sb

output = "./tests/pytest/output"
m = "~/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish"


def get_results(prefix, kb, oh):
    args = {
        "blast_results": f"{m}/1-First_pass/Unmatched/BLAST/unknown-blast.csv",
        "mapping": f"{m}/Databases/seq-header_mappings.tsv",
        "unknown_hits": f"{m}/1-First_pass/Combined/unknown_hits.tsv",
        "database_hits": f"{m}/1-First_pass/Combined/database_hits.tsv",
        "unmatched_fasta": f"{output}/{prefix}_unmatched.fasta",
        "unmatched_tsv": f"{output}/{prefix}_unmatched.tsv",
        "blast_query": "./tests/pytest/blast_query.txt",
        "unmatched_peptides": f"{m}/1-First_pass/Unmatched/unmatched_peptides.tsv",
        "output": f"{output}/{prefix}_blast_matched-test.tsv",
        "pep_threshold": 0.05,
        "evalue_threshold": 0.05,
        "identity_threshold": 0.99,
        "keep_best": kb,
        "one_hit": oh,
    }
    f = sb.main(args)
    return f


f = get_results("nO_nD", False, True)
assert len(f["ProteinId"]) == len(f["ProteinId"].unique())


# get_results("O_D", True, False)
# O_D_unmatched = pd.read_csv("./output/O_D_unmatched.tsv", sep="\t")
# #
# get_results("nO_D", False, False)
# nO_D_unmatched = pd.read_csv("./output/nO_D_unmatched.tsv", sep="\t")

# get_results("O_nD", False, True)
# O_nD_unmatched = pd.read_csv("./output/O_nD_unmatched.tsv", sep="\t")
