#!/usr/bin/env ipython

from Bio import SeqIO
from pyteomics import mass
import sys
import pandas as pd

db = sys.argv[1]
out = sys.argv[2]

fasta_frame = {"id": [], "seq": []}
for record in SeqIO.parse(db, format="fasta"):
    fasta_frame["id"].append(record.id)
    fasta_frame["seq"].append(str(record.seq))
fasta_frame = pd.DataFrame(fasta_frame)

def mass_calc(seq) -> float:
    if (intersect := {"X","B", "Z"} & set(seq)):
        for i in intersect:
            seq = seq.replace(i, "")
    return mass.fast_mass(seq)

fasta_frame["mass"] = fasta_frame["seq"].apply(mass_calc)
fasta_frame["length"] = fasta_frame["seq"].apply(len)
fasta_frame.to_csv(out, sep="\t", index=None)
