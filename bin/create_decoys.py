#!/usr/bin/env python
import sys
import pandas as pd
from pyteomics import mass
from Bio import SeqIO


def mass_calc(seq) -> float:
    if (intersect := {"X", "B", "Z"} & set(seq)):
        for i in intersect:
            seq = seq.replace(i, "")
    return mass.fast_mass(seq)


def write_lines(file: str, line_list: list) -> None:
    with open(file, "a") as a:
        for line in line_list:
            a.write(line)


input = sys.argv[1]
mapping_name = sys.argv[2]
# input = "./all.fasta"
# mapping_name = "./headers.tsv"

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
    if ("-DENOVO" in record.id or "-TRANSCRIPTOME" in record.id):
        seq_id = f"O{num_other}"
        num_other += 1
        is_download = False
    else:
        seq_id = f"P{num_download}"
        num_download += 1
        is_download = True
    mapping["id"].append(seq_id)
    mapping["header"].append(f"{record.description}")
    if "J" in record.seq or "U" in record.seq:
        seq = str(record.seq).replace("J", "").replace("U", "")
    else:
        seq = str(record.seq)
    mapping["seq"].append(seq)
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
