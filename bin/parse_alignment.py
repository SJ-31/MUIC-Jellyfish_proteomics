#!/usr/bin/env python
import polars as pl

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


def count_mismatch(
    alignment, id, protein_id, seq_mapping, result_dict, start, end
) -> None:
    cur_seq = seq_mapping[protein_id][start - 1 : end]
    cur_align = alignment[start - 1 : end]
    for index, chars in enumerate(zip(cur_seq, cur_align)):
        old, new = chars
        if old != new and new != "-":
            result_dict["ProteinId"].append(protein_id)
            result_dict["change"].append(f"{old}->{new}")
            result_dict["type"].append(f"{AAs[old]}->{AAs[new]}")
            result_dict["id"].append(id)
            result_dict["index"].append(index + (start - 1))


def main(args):
    alignment = pl.read_csv(args["alignment_path"], separator="\t")
    seq_map = pl.read_csv(args["seq_header_map"], separator="\t").filter(
        pl.col("id").is_in(alignment["ProteinId"])
    )
    id2seq = dict(zip(seq_map["id"], seq_map["seq"]))
    results = {"ProteinId": [], "change": [], "type": [], "id": [], "index": []}
    for x in alignment.iter_rows(named=True):
        count_mismatch(
            x["alignment"],
            x["id"],
            x["ProteinId"],
            id2seq,
            results,
            x["start"],
            x["end"],
        )

    all_mismatches = pl.DataFrame(results)
    all_metrics = (
        (
            all_mismatches.group_by("ProteinId")
            .agg(pl.len())
            .rename({"len": "n_mismatch"})
        )
        .join(seq_map, left_on="ProteinId", right_on="id")
        .select(["ProteinId", "n_mismatch", "header"])
    )
    return all_mismatches, all_metrics


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--alignment_path")
    parser.add_argument("-s", "--seq_header_map")
    parser.add_argument("-r", "--mismatch_tsv")
    parser.add_argument("-m", "--metric_tsv")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    mismatches, metrics = main(args)
    mismatches.write_csv(args["mismatch_tsv"], separator="\t", null_value="NA")
    metrics.write_csv(args["metric_tsv"], separator="\t", null_value="NA")
