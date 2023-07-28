#!/usr/bin/env python
import sys
import pandas as pd
from Bio import SeqIO

mapping: dict = {"key": [], "val": []}
mapping_name = sys.argv[3]
combined: str = sys.argv[2]
for num, record in enumerate(SeqIO.parse(sys.argv[1], "fasta")):
    mapping["key"].append(f"P{num}")
    mapping["val"].append(record.id)
    if "J" in record.seq or "U" in record.seq:
        seq = str(record.seq).replace("J", "").replace("U", "")
    else:
        seq = str(record.seq)
    with open("combined.fasta", "a") as a:
        a.write(f">P{num}" + "\n")
        a.write(str(seq) + "\n")
    with open(combined, "a") as a:
        a.write(f">P{num}" + "\n")
        a.write(str(seq) + "\n")
        a.write(f">rev_P{num}" + "\n")
        a.write(str(seq[::-1]) + "\n")
mapping = pd.DataFrame(mapping)
mapping.to_csv(mapping_name, sep="\t", header=False, index=False)
