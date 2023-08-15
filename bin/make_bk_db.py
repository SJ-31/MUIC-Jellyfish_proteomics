#!/usr/bin/env python
import pandas as pd
import numpy as np
from random import randint
from argparse import ArgumentParser
from Bio import SeqIO
parser = ArgumentParser(
    prog="bk_decoys",
    description="""
    Filter a database using Percolator protein identifications, populating it with
    decoys following Bern & Kil (2011)'s method
    """
)
parser.add_argument("percolator_valid")
parser.add_argument("percolator_decoys")
parser.add_argument("database")
parser.add_argument("output_file")

args = parser.parse_args()
matches = args.percolator_valid
decoys = args.percolator_decoys
old_db = args.database
output = args.output_file
# matches = "../results/test_manifest/First_pass/Percolator/comet_percolator_proteins.tsv"
# old_db = "/home/sc31/mapping.csv"
# decoys = "../results/test_manifest/First_pass/Percolator/comet_percolator_decoy_proteins.tsv"
# output = "/home/sc31/test_bk.fasta"

hits = []
for fasta in [matches, decoys]:
    ids_list = list(pd.read_csv(fasta, sep="\t")["ProteinId"])
    for group in ids_list:
        hits.extend(group.split(","))

map = (pd.read_csv(old_db, sep="\t", index_col=0))
matched = map.loc[hits]
matched["id"] = matched.index
num_valid = sum(matched["id"].str[0:4] != "rev_")
num_decoys = matched.shape[0] - num_valid
valid_prot = matched[matched["id"].str[0:4] != "rev_"].reset_index(drop=True)
extra = {"id": [], "seq": []}

rounds = num_valid - num_decoys
if rounds > 0:
    for r in range(rounds):
        num = randint(0, valid_prot.shape[0]-1)
        pick = valid_prot.iloc[num]
        extra["id"].append(f'rev_{pick["id"][1:]}_addedN{randint(0, 1000)}')
        extra["seq"].append(pick["seq"][::-1])
extra = pd.DataFrame(extra)
matched = pd.concat([extra, matched], axis=0)
new_fasta = []
for row in matched.iterrows():
    new_fasta.append(f'>{row[1]["id"]}\n{row[1]["seq"]}\n')
with open(output, "w") as o:
    o.write(''.join(new_fasta))
