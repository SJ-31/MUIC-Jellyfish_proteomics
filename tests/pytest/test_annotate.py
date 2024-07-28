#!/usr/bin/env python
import sys

dir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
sys.path.append(f"{dir}/bin")
sys.path.append(f"{dir}/tests/pytest")
import annotate as an
import test_verify_output as vo
import pandas as pd

pth = f"{dir}/results/C_indra/1-First_pass/"
dirname = f"{dir}/results/Unmatched/Database-annotated"
nftest = f"{dir}/tests/nf-test-out/annotate"
test = f"{dir}/tests/results/Unmatched/Database-annotated"
pytest = f"{dir}/tests/pytest/output"


def test_anno():
    import os

    args = {
        "input": f"{pth}/Unmatched/BLAST/C_indra_blast_matched.tsv",
        "previous_saved": f"{dir}/data/reference/C_indra_previous_annotations.tsv",
        # "annotate_extra": True,
    }
    f = an.anno(args)
    t = pd.read_csv("needs_annotating.tsv", sep="\t")
    os.remove("needs_annotating.fasta")
    os.remove("needs_annotating.tsv")
    return {
        "annotated_extra": f,
        "needs_annotating_extra": t,
    }


def test_merge_eggnog():
    dirname = f"{dir}/tests/results/Unmatched/Database-annotated"
    args = {
        "more_anno": f"{pytest}/annotated_extra.tsv",
        "eggnog_tsv": f"{dirname}/annotate_eggnog_unmatched/eggnog_matched.tsv",
    }
    f = an.mergeAnnotatedEggnog(args)
    return f


def test_merge_interpro():
    dirname = f"{pth}/Unmatched/Database-annotated"
    args = {
        "more_anno": f"{dirname}/jellyfish_downloads_anno-2.tsv",
        "interpro_query": f"{dirname}/annotate_interpro_unmatched/needs_annotating2.fasta",
        "r_source": "./bin",
        "input": f"{dirname}/annotate_interpro_unmatched/annotated.tsv",
    }
    f = an.mergeAnnotatedInterpro(args)
    return f


a = test_anno()
# # Try drop duplicates by protein id
for name, df in a.items():
    print(f"{name}: {df.shape}")
    print(f"    Num unique protein ids {df['ProteinId'].unique().shape}")
    print(vo.columnsWith(",", df))
