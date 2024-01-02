#!/usr/bin/env python
# Obtain the sequence of eggnog proteins, writing them into the metadata file
import numpy as np
import pandas as pd
import subprocess


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--eggnog_database")
    parser.add_argument("-m", "--metadata_file")
    parser.add_argument("-o", "--output_file")
    parser.add_argument("-a", "--annotations_file")
    args = vars(parser.parse_args())  # convert to dict
    return args


def main(args: dict):
    query_file = "./queries.txt"
    anno = pd.read_csv(args["annotations_file"], sep="\t")
    eggnog_queries = "\n".join(list(anno["seed_ortholog"]))
    with open("./queries.txt", "w") as q:
        q.write(eggnog_queries)
    command = f'seqkit grep -n -r -f {query_file} {args["eggnog_database"]}'
    seqkit = subprocess.run(command, shell=True, capture_output=True)
    unformatted = str(seqkit.stdout)[2:-1].split("\\n")
    seqs = {"seed_ortholog": set(), "SO_seq": []}
    seq_count = 0
    for line in unformatted:
        if line.startswith(">"):
            seqs["seed_ortholog"].add(line[1:])
            seqs["SO_seq"].append("")
            seq_count += 1
        else:
            seqs["SO_seq"][seq_count - 1] = (
                seqs["SO_seq"][seq_count - 1] + line
            )
    for so in anno["seed_ortholog"]:
        if so not in seqs["seed_ortholog"]:
            seqs["seed_ortholog"].add(so)
            seqs["SO_seq"].append(np.NaN)
    seqs["seed_ortholog"] = list(seqs["seed_ortholog"])
    seqs = pd.DataFrame(seqs)
    seqs.to_csv("sequence_df.tsv", sep="\t", index=False)
    anno = anno.filter(["ProteinId", "seed_ortholog"]).merge(
        seqs, on="seed_ortholog"
    )
    meta = pd.read_csv(args["metadata_file"], sep="\t")
    final = anno.filter(["ProteinId", "SO_seq"]).merge(meta, on="ProteinId")
    final.to_csv(args["output_file"], sep="\t", na_rep="NA", index=None)


if __name__ == "__main__":
    args = parse_args()
    main(args)
