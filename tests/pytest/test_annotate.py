#!/usr/bin/env ipython
import sys

sys.path.append("./bin")
import annotate as an
import pandas as pd

pth = "./results/jellyfish/1-First_pass/"
dirname = "./results/Unmatched/Database-annotated"
args = {
    "r_source": "./bin",
    "eggnog_anno_tsv": f"{dirname}/anno-eggnog.tsv",
    "eggnog_meta_tsv": f"{dirname}/meta-eggnog.tsv",
    "meta_tsv": f"{dirname}/O_nD_meta.tsv",
    "anno_tsv": f"{dirname}/O_nD_anno.tsv",
    "interpro_query": f"{dirname}/needs_annotating2.fasta",
    "input": f"{pth}/Unmatched/BLAST/db_hits-O_nD.tsv",
}

f = an.anno(args)
