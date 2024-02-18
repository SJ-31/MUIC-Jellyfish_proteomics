#!/usr/bin/env python

import re
from Bio import SeqIO


def cleanPeptide(peptide):
    peptide = peptide.upper().replace("X", "")
    return "".join(re.findall("[A-Z]+", peptide))


def cleanFasta(file):
    cleaned_fasta = ""
    for record in SeqIO.parse(file, "fasta"):
        if "[" in record.seq:
            string = f">{record.id}\n{cleanPeptide(str(record.seq))}\n"
        else:
            string = f">{record.id}\n{str(record.seq)}\n"
        cleaned_fasta = cleaned_fasta + string
    return cleaned_fasta


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    cleaned = cleanFasta(args["input"])
    with open(args["output"], "w") as w:
        w.write(cleaned)
