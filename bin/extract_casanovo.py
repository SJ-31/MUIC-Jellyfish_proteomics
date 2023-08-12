#!/usr/bin/env ipython
import pandas as pd
import sys
import re

file = pd.read_csv(sys.argv[1], sep="\t")
peps = file["sequence"]

def get_seq(index: int, peptide: str) -> list[str]:
    peps: list[str] = [pep for pep in re.findall(r"([A-Z]+)", peptide)]
    rev = ''.join(peps)[::-1]
    return [f">D{index}\n"] + peps + ["\n"] + [f">rev_D{index}\n{rev}\n"]

w_decoys: list = []
decoys: list = []
normal: list = []
for index, p in enumerate(peps):
    seq = ''.join(re.findall(r"[A-Z]+", p))
    pep = f">D{index}\n{seq}\n"
    rev = f">rev_D{index}\n{seq[::-1]}\n"
    w_decoys.extend([pep, rev])
    normal.append(pep)
    decoys.append(rev)

output = ["casanovo_combined.fasta", "casanovo_decoys.fasta", "casanovo_normal.fasta"]
for file_name, out in zip(output, [w_decoys, decoys, normal]):
    with open(file_name, "w") as f:
        f.write(''.join(out))
