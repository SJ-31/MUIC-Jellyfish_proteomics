#!/usr/bin/env python
import random
import sys
import pandas as pd
from pyteomics import mass
from Bio import SeqIO


def resolveResidue(seq):
    non_standard: dict = {
        "X": [""],
        "B": ["N", "D"],
        "Z": ["Q", "E"],
        "J": ["L", "I"],
        "U": [""],
    }
    if intersect := set(non_standard.keys()) & set(seq):
        for i in intersect:
            seq = seq.replace(i, random.choice(non_standard[i]))
    return seq


def mass_calc(seq) -> float:
    return mass.fast_mass(resolveResidue(seq))


def write_lines(file: str, line_list: list) -> None:
    with open(file, "a") as a:
        for line in line_list:
            a.write(line)


input = sys.argv[1]
mapping_name = sys.argv[2]

mapping: dict = {
    "id": [],
    "header": [],
    "seq": [],
}
combined: str = "decoysWnormal.fasta"
decoys: str = "all_decoys.fasta"
downloaded: str = "downloaded.fasta"
normal: str = "all_normal.fasta"
num_download = 0
num_other = 0
is_download = True
for record in SeqIO.parse(input, "fasta"):
    if "-DENOVO" in record.id:
        seq_id = f"D{num_other}"
        num_other += 1
        is_download = False
    elif "-TRANSCRIPTOME" in record.id:
        seq_id = f"T{num_other}"
        num_other += 1
        is_download = False
    else:
        seq_id = f"P{num_download}"
        num_download += 1
        is_download = True
    mapping["id"].extend([seq_id, f"rev_{seq_id}"])
    mapping["header"].extend([f"{record.description}", "DECOY"])
    seq = resolveResidue(str(record.seq))
    mapping["seq"].extend([seq, str(seq[::-1])])
    norm_lines: list = [f">{seq_id}" + "\n", str(seq) + "\n"]
    decoy_lines: list = [f">rev_{seq_id}" + "\n", str(seq[::-1]) + "\n"]
    write_lines(normal, norm_lines)
    write_lines(decoys, decoy_lines)
    write_lines(combined, norm_lines + decoy_lines)
    if is_download:
        write_lines(downloaded, norm_lines)
        is_download = False

mapping = pd.DataFrame(mapping)
mapping["mass"] = mapping["seq"].apply(mass_calc)
mapping["length"] = mapping["seq"].apply(len)
mapping.to_csv(mapping_name, sep="\t", index=False)
