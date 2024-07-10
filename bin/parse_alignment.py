#!/usr/bin/env python
import pandas as pd
from Bio import SeqIO

AAs = {
    "A": "NP",
    "G": "NP",
    "L": "NP",
    "I": "NP",
    "W": "NP",
    "F": "NP",
    "V": "NP",
    "P": "NP",
    "M": "NP",
    "H": "PB",
    "K": "PB",
    "R": "PB",
    "E": "PA",
    "D": "PA",
    "Y": "PN",
    "C": "PN",
    "T": "PN",
    "U": "PA",
    "N": "PN",
    "O": "PB",
    "Q": "PN",
    "S": "PN",
}


def record_mismatch(
    protein_id: str, old: str, new: str, index: int
) -> pd.DataFrame:
    temp = {"ProteinId": [], "change": [], "type": [], "index": []}
    temp["ProteinId"].append(protein_id)
    temp["change"].append(f"{old}->{new}")
    temp["type"].append(f"{AAs[old]}->{AAs[new]}")
    temp["index"].append(index)
    row = pd.DataFrame(temp)
    return row


def record_alignment(header: str, seq, alignment) -> tuple:
    protein_id = header[: header.find(":")]
    header = header[header.find(":") + 1 :]
    if "DENOVO" in header or "U" in protein_id:
        return tuple()
    length = len(seq)
    matches: int = length
    n_mismatches = 0
    i: int = 0
    mismatches = pd.DataFrame()
    while i < length:
        old: str = seq[i]
        new: str = alignment[i]
        if new == "[":
            indices = (i + 1, alignment[i:].find("]") + i)
            possible_mismatches = alignment[indices[0] : indices[1]]
            n_mismatches += 1
            for r in possible_mismatches:
                if old != r:
                    mismatches = pd.concat(
                        [
                            mismatches,
                            record_mismatch(protein_id, old, r, i),
                        ]
                    )
            alignment = (
                alignment[: indices[0] - 1] + "_" + alignment[indices[1] + 1 :]
            )
        elif new == "-":
            matches -= 1
        elif old != new:
            mismatches = pd.concat(
                [mismatches, record_mismatch(protein_id, old, new, i)]
            )
            n_mismatches += 1
        i += 1
    metrics = {
        "ProteinId": protein_id,
        "header": header,
        "n_mismatches": n_mismatches,
        "pcoverage_nmatch": matches / length,
    }
    return mismatches, pd.DataFrame(metrics, index=[0])


def parse_alignment_file(file) -> tuple:
    lines = SeqIO.parse(file, "fasta")
    metric_df = pd.DataFrame()
    mismatch_df = pd.DataFrame()
    for l in lines:
        header, seq = l.id, l.seq
        if "ALIGNED" not in header:
            alignment = next(lines)
            result = record_alignment(header, seq, alignment.seq)
            if not result:
                continue
            mismatch_df = pd.concat([mismatch_df, result[0]])
            metric_df = pd.concat([metric_df, result[1]])
    return metric_df, mismatch_df


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-r", "--mismatch_tsv")
    parser.add_argument("-m", "--metric_tsv")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    m, r = parse_alignment_file(args["input"])
    r.to_csv(args["mismatch_tsv"], index=False, sep="\t")
    m.to_csv(args["metric_tsv"], index=False, sep="\t")
