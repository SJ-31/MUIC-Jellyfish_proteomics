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
parser.add_argument("mapping_file")
parser.add_argument("output_file")

args = parser.parse_args()
matches = args.percolator_valid
decoys = args.percolator_decoys
old_db = args.database
output = args.output_file
map = args.mapping_file
# matches = "/home/sc31/Bio_SDD/MUIC_senior_project/workflow/results/test_manifest/testing_intersect/comet_percolator_proteins.tsv"
# map = "/home/sc31/Bio_SDD/MUIC_senior_project/workflow/results/test_manifest/testing_intersect/header_mappings.tsv"
# old_db = "/home/sc31/Bio_SDD/MUIC_senior_project/workflow/results/test_manifest/testing_intersect/decoysWnormal.fasta"
# decoys = "/home/sc31/Bio_SDD/MUIC_senior_project/workflow/results/test_manifest/Percolator/comet_percolator_decoy_proteins.tsv"
# output = "/home/sc31/test_bk.fasta"

hits = []
for fasta in [matches, decoys]:
    hits.extend(list(pd.read_csv(fasta, sep="\t")["ProteinId"]))
hits = np.array(hits)
num_valid: int = 0
new_db = {"id": [], "seq": [], "is_decoy": []}

for record in SeqIO.parse(old_db, format="fasta"):
    if (id := record.id) in hits:
        new_db["id"].append(f">{id}")
        new_db["seq"].append(str(record.seq))
        if record.id[0:4] != "rev_":
            new_db["is_decoy"].append(False)
            num_valid += 1
        else:
            new_db["is_decoy"].append(True)
    if len(new_db["id"]) == len(hits):
        break
new_db = pd.DataFrame(new_db)
valid_prot = new_db[new_db["is_decoy"] == False].reset_index()
extra = {"id": [], "seq": [], "is_decoy": []}

rounds = num_valid - (new_db.shape[0] - num_valid)
print(rounds)
if rounds > 0:
    for r in range(rounds):
        num = randint(0, valid_prot.shape[0]-1)
        pick = valid_prot.iloc[num]
        extra["id"].append(f'>rev_{pick["id"][1:]}_addedN{randint(0, 1000)}')
        extra["seq"].append(pick["seq"][::-1])
        extra["is_decoy"].append(True)
extra = pd.DataFrame(extra)
new_db = pd.concat([extra, new_db], axis=0)
new_fasta = []
for row in new_db.iterrows():
    new_fasta.append(f'{row[1]["id"]}\n{row[1]["seq"]}\n')
with open(output, "w") as o:
    o.write(''.join(new_fasta))
