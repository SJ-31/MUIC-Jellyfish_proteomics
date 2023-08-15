#!/usr/bin/env ipython

from Bio import SeqIO
import sys
import pandas as pd

db = sys.argv[1]
out = sys.argv[2]

fasta_frame = {"id": [], "seq": []}
for record in SeqIO.parse(db, format="fasta"):
    fasta_frame["id"].append(record.id)
    fasta_frame["seq"].append(record.seq)
fasta_frame = pd.DataFrame(fasta_frame)
fasta_frame.to_csv(out, sep="\t", index=None)
