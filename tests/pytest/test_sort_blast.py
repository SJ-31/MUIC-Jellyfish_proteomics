#!/usr/bin/env python
import os
import sys


workdir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")

import sort_blast as sb

output = f"{workdir}/tests/pytest/output"
m = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish"


def get_results():
    args = {
        "blast_results": f"{m}/1-First_pass/Unmatched/BLAST/unknown-blast.csv",
        "mapping": f"{m}/Databases/seq-header_mappings.tsv",
        "unknown_hits": f"{m}/1-First_pass/Combined/unknown_hits.tsv",
        "database_hits": f"{m}/1-First_pass/Combined/database_hits.tsv",
        "unmatched_fasta": f"{m}/1-First_pass/Unmatched/BLAST/jellyfish_blast_unmatched.fasta",
        "unmatched_tsv": f"{m}/1-First_pass/Unmatched/BLAST/jellyfish_blast_unmatched.tsv",
        "blast_query": f"{workdir}/tests/results/Unmatched/queries.txt",
        "unmatched_peptides": f"{m}/1-First_pass/Unmatched/unmatched_peptides.tsv",
        "adjust": True,
        "pep_threshold": 0.05,
        "evalue_threshold": 0.05,
        "identity_threshold": 0.99,
    }
    f = sb.main(args)
    return f


f = get_results()
assert len(f["ProteinId"]) == len(f["ProteinId"].unique())
assert columnsWith(",", f) == ["header"]
has_commas = f["peptideIds"][f["peptideIds"].str.contains(",")]
unknown_hits = pd.read_csv(
    f"{m}/1-First_pass/Combined/unknown_hits.tsv", sep="\t"
)
db_hits = pd.read_csv(f"{m}/1-First_pass/Combined/database_hits.tsv", sep="\t")
peps = pd.read_csv(
    f"{m}/1-First_pass/Unmatched/unmatched_peptides.tsv", sep="\t"
)
columnsWith(",", peps)
columnsWith(",", unknown_hits)
columnsWith(",", db_hits)
