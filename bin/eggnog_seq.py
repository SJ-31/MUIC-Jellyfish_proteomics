#!/usr/bin/env python
# Obtain the sequence of eggnog proteins, writing them into the metadata file
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
    eggnog_queries = '\n'.join(list(anno["seed_ortholog"]))
    with open("./queries.txt", "w") as q:
        q.write(eggnog_queries)
    command = f'seqkit grep -n -r -f {query_file} {args["eggnog_database"]}'
    seqkit = subprocess.run(command, shell=True, capture_output=True)
    unformatted = str(seqkit.stdout)[2:-1].split("\\n")
    seqs = {"seed_ortholog": [], "seq": []}
    seq_count = 0
    for line in unformatted:
        if line.startswith(">"):
            seqs["seed_ortholog"].append(line[1:])
            seqs["seq"].append("")
            seq_count += 1
        else:
            seqs["seq"][seq_count-1] = (seqs["seq"][seq_count-1] + line)
    seqs = pd.DataFrame(seqs)
    anno = (anno.filter(["ProteinId", "seed_ortholog"])
            .merge(seqs, on="seed_ortholog"))
    meta = pd.read_csv(args["metadata_file"], sep="\t").drop("seq",
                                                             axis="columns")
    final = anno.filter(["ProteinId", "seq"]).merge(meta, on="ProteinId")
    final.to_csv(args["output_file"], sep="\t", na_rep="-", index=None)


if __name__ == '__main__':
    args = parse_args()
    main(args)
