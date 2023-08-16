#!/usr/bin/env python
import sys
import pandas as pd
from Bio import SeqIO

def write_lines(file: str, line_list: list) -> None:
    with open(file, "a") as a:
        for line in line_list:
            a.write(line)

mapping: dict = {"key": [], "val": []}
mapping_name = sys.argv[2]
combined: str = "decoysWnormal.fasta"
decoys: str = "all_decoys.fasta"
normal: str = "all_normal.fasta"
for num, record in enumerate(SeqIO.parse(sys.argv[1], "fasta")):
    mapping["key"].append(f"P{num}")
    mapping["val"].append(f"{record.description}")
    if "J" in record.seq or "U" in record.seq:
        seq = str(record.seq).replace("J", "").replace("U", "")
    else:
        seq = str(record.seq)
    norm_lines: list = [f">P{num}" + "\n", str(seq) + "\n"]
    decoy_lines: list = [f">rev_P{num}" + "\n", str(seq[::-1]) + "\n"]
    write_lines(normal, norm_lines)
    write_lines(decoys, decoy_lines)
    write_lines(combined, norm_lines + decoy_lines)
mapping = pd.DataFrame(mapping)
mapping.to_csv(mapping_name, sep="\t", header=False, index=False)
