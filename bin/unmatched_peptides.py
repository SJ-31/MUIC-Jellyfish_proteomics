#!/usr/bin/env python
"""
Obtain peptides that have not been matched to proteins in the database
search/protein inference with percolator, but pass a q-value and pep
threshold
"""

from pathlib import Path

import re
import pandas as pd
import numpy as np


def perc_row(series_row):
    return pd.DataFrame(
        {
            "PSMId": series_row[0],
            "score": series_row[1],
            "q-value": series_row[2],
            "posterior_error_prob": series_row[3],
            "peptide": series_row[4],
            "proteinIds": ";".join(series_row[5:]),
        },
        index=[1],
    )


def read_percolator(filepath, p_threshold, q_threshold):
    percolator = (
        pd.read_csv(filepath, index_col=False)
        .iloc[:, 0]
        .apply(str.split, sep="\t")
        .apply(perc_row)
        .to_list()
    )
    final = pd.concat(percolator)
    final["posterior_error_prob"] = final["posterior_error_prob"].astype(float)
    final["q-value"] = final["q-value"].astype(float)
    final = final[final["q-value"] <= q_threshold]
    final = final[final["posterior_error_prob"] <= p_threshold]
    final["peptide"] = final["peptide"].apply(
        lambda x: x.upper()[2:-2]
    )  # Remove N- and C-termini, formatting peptide string like in the
    # percolator "proteins" file
    return final


def read_tide(filepath, p_threshold, q_threshold):
    tide = pd.read_csv(filepath, sep="\t")
    tide["sequence"] = tide["sequence"]
    tide = tide[tide["percolator q-value"] <= q_threshold]
    tide = tide[tide["percolator PEP"] <= p_threshold].rename(
        columns={"sequence": "peptide", "protein id": "proteinIds"}
    )
    return tide


def get_engine_files(path) -> dict:
    files = [file.absolute() for file in Path(path).glob("*_percolator_*")]
    names = [re.sub(r"_percolator_.*", "", f.name) for f in files]
    engine_files = dict(zip(names, files))
    return engine_files


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_path")
    parser.add_argument("-t", "--unmatched_tsv")
    parser.add_argument("-p", "--pep_threshold")
    parser.add_argument("-q", "--q_threshold")
    args = vars(parser.parse_args())  # convert to dict
    return args


def main(args: dict):
    all_unmatched = set()
    unmatched_df = []
    pep_threshold = float(args["pep_threshold"])
    q_threshold = float(args["q_threshold"])
    path = args["input_path"]
    engine_files = get_engine_files(path)
    for engine, path in engine_files.items():
        if engine == "tide":
            current = read_tide(
                engine_files["tide"], pep_threshold, q_threshold
            )
        else:
            current = read_percolator(path, pep_threshold, q_threshold)
        unmatched = pd.Series(
            current.where(current["proteinIds"] == "")["peptide"].unique()
        )
        unmatched_df.append(current[current["peptide"].isin(unmatched)])
        all_unmatched = all_unmatched | set(unmatched)
    all_unmatched = pd.Series(list(all_unmatched)).reset_index()
    all_unmatched["ProteinId"] = all_unmatched["index"].apply(
        lambda x: f"U{x}"
    )
    all_unmatched.drop("index", axis="columns", inplace=True)
    all_unmatched.rename({0: "seq"}, axis="columns", inplace=True)
    all_unmatched_df = pd.concat(unmatched_df).groupby("peptide").sample()
    all_unmatched_df.rename(
        {"q-value": "q.value", "peptide": "peptideIds"},
        axis="columns",
        inplace=True,
    )
    df = all_unmatched.merge(
        all_unmatched_df, how="inner", left_on="seq", right_on="peptideIds"
    )
    df["ProteinGroupId"] = np.NaN
    df["Group"] = np.NaN
    df["header"] = "unknown"
    df = df[
        [
            "ProteinId",
            "ProteinGroupId",
            "q.value",
            "posterior_error_prob",
            "peptideIds",
            "header",
            "Group",
        ]
    ]
    df.to_csv(args["unmatched_tsv"], sep="\t", index=False, na_rep="-")


if __name__ == "__main__":
    args = parse_args()
    main(args)
