#!/usr/bin/env ipython

import pandas as pd
test = "comet_percolator_proteins.tsv"
table = pd.read_csv(test, sep="\t")
sequence =


def remove_decoys(prot_list):
    prots = prot_list.split(",")
    return ','.join(p for p in prots if "_" not in p)

table["ProteinId"] = table["ProteinId"].apply(remove_decoys)
