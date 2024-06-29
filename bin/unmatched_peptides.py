#!/usr/bin/env python
"""
Obtain peptides that have not been matched to proteins in the database
search/protein inference with percolator, but pass a q-value and pep
threshold
"""

from pathlib import Path

import re
import pandas as pd
import polars as pl


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
    files = [file.absolute() for file in Path(path).glob("*_percolator_psms.tsv")]
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
    unmatched_list = []
    pep_threshold = float(args["pep_threshold"])
    q_threshold = float(args["q_threshold"])
    path = args["input_path"]
    engine_files = get_engine_files(path)
    for engine, path in engine_files.items():
        if engine == "tide":
            current = read_tide(engine_files["tide"], pep_threshold, q_threshold)
        else:
            current = read_percolator(path, pep_threshold, q_threshold)
        cur_df = current.where(current["proteinIds"] == "").drop_duplicates("peptide")
        cur_df["engine"] = engine
        unmatched_list.append(cur_df)

    unmatched_df = pl.from_pandas(pd.concat(unmatched_list))
    unmatched_df = (
        unmatched_df.group_by("peptide")
        .agg(
            pl.col("q-value").mean(),
            pl.col("posterior_error_prob").mean(),
            pl.col("engine"),
        )
        .with_columns(pl.col("engine").list.join(";"))
    )
    new_ids = pl.Series([f"U{i}" for i in range(unmatched_df.shape[0])])
    unmatched_df = (
        unmatched_df.with_columns(
            ProteinId=new_ids,
            header=pl.lit("unknown"),
            ProteinGroupId=pl.lit("U"),
            Group=pl.lit("U"),
        )
        .rename({"q-value": "q.value", "peptide": "peptideIds"})
        .select(
            [
                "ProteinId",
                "ProteinGroupId",
                "q.value",
                "posterior_error_prob",
                "peptideIds",
                "header",
                "engine",
                "Group",
            ]
        )
    )

    return unmatched_df


if __name__ == "__main__":
    args = parse_args()
    result = main(args)
    result.write_csv(args["unmatched_tsv"], separator="\t")
