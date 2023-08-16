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

normal: list = []
for index, p in enumerate(peps):
    seq = ''.join(re.findall(r"[A-Z]+", p))
    pep = f">{name}{index}-DENOVO\n{seq}\n"
    normal.append(pep)

with open(f"{name.lower()}_normal.fasta", "w") as f:
    f.write(''.join(normal))
