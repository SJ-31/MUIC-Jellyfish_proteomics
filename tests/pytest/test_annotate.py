#!/usr/bin/env ipython
import sys

sys.path.append("./bin")
import annotate as an
import pandas as pd

pth = "./results/jellyfish/1-First_pass/"
dirname = "./results/Unmatched/Database-annotated"
nftest = "./tests/nf-test-out/annotate"
test = "./tests/results/Unmatched/Database-annotated"
pytest = "./tests/pytest/output"


def test_anno():
    args = {
        "input": f"{pth}/Unmatched/BLAST/jellyfish_blast_matched.tsv",
    }
    f = an.anno(args)
    return f


def test_merge_eggnog():
    dirname = "./tests/results/Unmatched/Database-annotated"
    args = {
        "more_anno": f"{dirname}/jellyfish_downloads_anno-1.tsv",
        "eggnog_tsv": f"{dirname}/annotate_eggnog_unmatched/eggnog_matched.tsv",
    }
    f = an.mergeAnnotatedEggnog(args)
    return f


def test_merge_interpro():
    dirname = "./tests/results/Unmatched/Database-annotated"
    args = {
        "more_anno": f"{dirname}/jellyfish_downloads_anno-2.tsv",
        "interpro_query": f"{dirname}/annotate_interpro_unmatched/needs_annotating2.fasta",
        "r_source": "./bin",
        "input": f"{dirname}/annotate_interpro_unmatched/annotated.tsv",
    }
    f = an.mergeAnnotatedInterpro(args)
    return f


f = test_anno()
f.to_csv(
    f"{pytest}/jellyfish_downloads_anno-1.tsv",
    sep="\t",
    index=False,
    na_rep="NaN",
)
eggnog = test_merge_eggnog()
i = test_merge_interpro()
