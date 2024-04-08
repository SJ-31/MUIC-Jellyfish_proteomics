#!/usr/bin/env python

import pandas as pd


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-s", "--sequence_col")
    parser.add_argument("-d", "--header_col")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parseArgs()
    if ".tsv" in args["input"]:
        df = pd.read_csv(args["input"], sep="\t")
    else:
        df = pd.read_csv(args["input"])
    seqs = df[args["sequence_col"]]
    headers = df[args["header_col"]]
    fasta = []
    for h, s in zip(headers, seqs):
        fasta.append(f">{h}\n{s}")
    fasta.append("\n")
    with open(args["output"], "w", encoding="utf-8") as f:
        f.write("\n".join(fasta))
