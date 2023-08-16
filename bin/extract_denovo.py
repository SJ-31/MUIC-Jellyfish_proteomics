#!/usr/bin/env ipython
import pandas as pd
import sys
import re

file = pd.read_csv(sys.argv[1], sep="\t")
engine = sys.argv[2]
if engine == "casanovo":
    peps = file["sequence"]
    name = "Casanovo"
elif engine == "pepnet":
    peps = file["DENOVO"]
    name = "PepNet"


w_decoys: list = []
decoys: list = []
normal: list = []
for index, p in enumerate(peps):
    seq = ''.join(re.findall(r"[A-Z]+", p))
    pep = f">{name}{index}-DENOVO\n{seq}\n"
    rev = f">rev_{name}{index}-DENOVO\n{seq[::-1]}\n"
    w_decoys.extend([pep, rev])
    normal.append(pep)
    decoys.append(rev)

output = [f"{name.lower()}_combined.fasta", f"{name.lower()}_decoys.fasta", f"{name.lower()}_normal.fasta"]
for file_name, out in zip(output, [w_decoys, decoys, normal]):
    with open(file_name, "w") as f:
        f.write(''.join(out))
