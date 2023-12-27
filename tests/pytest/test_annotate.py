#!/usr/bin/env ipython
import pandas as pd

# interpro = ./Unmatched/InterPro/no_one_hits_degenerates_unmatched_eggnog-SCAN.tsv
# sim = unannotated["ProteinId"].sample(sorted[0:54].shape[0])
dirname = "../results/Unmatched/Database-annotated"
args = {
    "r_source": "../../bin",
    "eggnog_anno_tsv": f"{dirname}/anno-eggnog.tsv",
    "eggnog_meta_tsv": f"{dirname}/meta-eggnog.tsv",
    "meta_tsv": f"{dirname}/O_nD_meta.tsv",
    "anno_tsv": f"{dirname}/O_nD_anno.tsv",
    "interpro_query": f"{dirname}/needs_annotating2.fasta",
    "input": f"{dirname}/db_hits-O_nD.tsv",
}

unmatched = f"{dirname}/needs_annotating.fasta"
with open(unmatched, "r") as u:
    headers = [u.strip().replace(">", "") for u in u.readlines() if ">" in u]
    headers = pd.Series(headers)

meta = pd.read_csv(args["meta_tsv"], sep="\t")
anno = pd.read_csv(args["anno_tsv"], sep="\t")
final = anno.merge(meta, on="ProteinId")
final[final["ProteinId"].isin(headers)].to_csv(
    "needs_annotating.tsv", sep="\t", index=False
)
