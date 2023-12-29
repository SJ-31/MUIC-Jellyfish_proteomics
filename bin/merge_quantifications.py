#!/usr/bin/env python

import numpy as np
import pandas as pd


def clean(target, string):
    target = set(target)
    if string in target:
        target.remove(string)
    return list(target)


def read_directlfq(directlfq):
    df = pd.read_csv(directlfq, sep="\t").iloc[:, 1:]
    rename_mapping = dict(
        zip(df.columns[1:], [f"directlfq-{col}" for col in df.columns[1:]])
    )
    df = df.rename(columns=rename_mapping)
    return df


def clean_list_col(colname, frame, split):
    return (
        frame[colname]
        .apply(str.split, sep=split)
        .apply(clean, string="")
        .apply(clean, string="NA")
    )


def sort_val(vals):
    if len(vals) > 1:
        # Sort the VALs and return the second smallest val
        # If BOTH engines have small vals for the protein,
        #   then it will be accepted
        sorted_vals = [float(vals) for vals in vals]
        sorted_vals.sort()
        return sorted_vals[1]
    return float(vals[0])


def split_proteins(row):
    new_len = len(row["protein"])
    new = row.drop("protein").map(lambda x: [x] * new_len)
    new["protein"] = row["protein"]
    new_df = pd.DataFrame(new.to_dict())
    new_df = new_df[[new_df.columns[-1]] + new_df.columns[:-1].to_list()]
    return new_df


def merge_dlfq(direct_lfq, prot_df, pep_thresh, fdr):
    dlfq_copy = direct_lfq.copy()
    dlfq_copy["protein"] = clean_list_col("protein", dlfq_copy, ";").apply(
        list
    )
    dlfq_copy = dlfq_copy[dlfq_copy["protein"].map(lambda x: len(x) > 0)]
    split = [split_proteins(row[1]) for row in dlfq_copy.iterrows()]
    dlfq_copy = pd.concat(split)
    dlfq_copy = dlfq_copy[~dlfq_copy["protein"].str.contains("rev_")]
    dlfq_copy = dlfq_copy.groupby("protein").median()
    merged = pd.merge(
        prot_df, dlfq_copy, how="left", left_on="ProteinId", right_on="protein"
    )
    merged["posterior_error_prob"] = clean_list_col(
        "posterior_error_prob", merged, ","
    ).apply(sort_val)
    merged["q.value"] = clean_list_col("q.value", merged, ",").apply(sort_val)
    merged = merged[merged["posterior_error_prob"] < pep_thresh]
    merged = merged[merged["q.value"] < fdr]
    return merged


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--direct_lfq")
    parser.add_argument("-f", "--flash_lfq")
    parser.add_argument("-i", "--intersected_searches")
    parser.add_argument("-s", "--open_searches")
    parser.add_argument("-u", "--unmatched_peptides")
    parser.add_argument("-p", "--pep_threshold")
    parser.add_argument("-q", "--fdr_threshold")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


def concatComma(df, col1, col2) -> pd.Series:
    """
    Concatenate columns col1 and col2 such that their entries are
    joined end-to-end in a csv list.
    """
    return df[col1].combine(df[col2], lambda x, y: f"{x},{y}")


def read_flashlfq(file):
    protdf = (
        pd.read_csv(file, sep="\t")
        .drop(0)
        .drop(["Gene Name", "Organism"], axis="columns")
        .query('~(`Protein Groups`.str.startswith("rev_"))')
    )
    intensity_cols = protdf.loc[
        :, protdf.columns != "Protein Groups"
    ].columns.to_series()
    new_cols = list(
        intensity_cols.apply(
            lambda x: f"flashlfq_{x.replace('Intensity_Default_', 'intensity')}"
        )
    )
    name_mapper = {old: new for new, old in zip(new_cols, intensity_cols)}
    protdf.rename(name_mapper, axis="columns", inplace=True)
    return protdf


def main(args: dict):
    dlfq = read_directlfq(args["direct_lfq"])
    flfq = read_flashlfq(args["flash_lfq"])
    proteins = pd.read_csv(args["intersected_searches"], sep="\t")
    open_searches = pd.read_csv(args["open_searches"], sep="\t")
    open_searches.rename({"q-value": "q.value"}, axis="columns", inplace=True)

    # Identify proteins that were found in both standard and open search
    shared = proteins.merge(
        open_searches.filter(
            ["ProteinId", "peptideIds", "q.value", "posterior_error_prob"]
        ),
        on="ProteinId",
        how="inner",
    )
    shared["peptideIds"] = concatComma(shared, "peptideIds_x", "peptideIds_y")
    shared["q.value"] = concatComma(shared, "q.value_x", "q.value_y")
    shared["posterior_error_prob"] = concatComma(
        shared, "posterior_error_prob_x", "posterior_error_prob_y"
    )
    shared.drop(
        [
            "peptideIds_x",
            "peptideIds_y",
            "posterior_error_prob_x",
            "posterior_error_prob_y",
            "q.value_x",
            "q.value_y",
        ],
        inplace=True,
        axis="columns",
    )

    # Isolate proteins found only in open search and standard search
    proteins = proteins[~proteins["ProteinId"].isin(shared["ProteinId"])]
    open_searches = open_searches[
        ~open_searches["ProteinId"].isin(shared["ProteinId"])
    ]

    # Merge with quantification
    merged = merge_dlfq(
        dlfq,
        proteins,
        float(args["pep_threshold"]),
        float(args["fdr_threshold"]),
    )
    merged_shared = merge_dlfq(
        dlfq,
        shared,
        float(args["pep_threshold"]),
        float(args["fdr_threshold"]),
    )
    merged_open = merge_dlfq(
        dlfq,
        open_searches,
        float(args["pep_threshold"]),
        float(args["fdr_threshold"]),
    )
    final = pd.concat([merged, merged_shared, merged_open])
    final["Anno_method"] = "initial_database"
    final = final.merge(
        flfq, left_on="ProteinId", right_on="Protein Groups", how="left"
    ).drop("Protein Groups", axis="columns")
    return final


def check_unmatched(final_frame, peps):
    """
    Filter out the peptides in the dataframe "peps" to leave only peptides that
    are truly unmatched
    """
    pep_df = pd.read_csv(peps, sep="\t")
    all_peps = pd.Series(
        [
            pep
            for pep_list in final_frame["peptideIds"].str.split(",")
            for pep in pep_list
        ]
    )
    actual_unmatched = pep_df.query("~(`peptideIds`.isin(@all_peps))")
    fasta_str = "".join(
        actual_unmatched.apply(
            lambda x: f'>{x["ProteinId"]}\n{x["peptideIds"]}\n', axis="columns"
        ).to_list()
    )
    return (actual_unmatched, fasta_str)


if __name__ == "__main__":
    args = parse_args()
    m = main(args)
    m.to_csv(args["output"], sep="\t", index=False)
    unmatched_peptides = check_unmatched(m, args["unmatched_peptides"])
    unmatched_peptides[0].to_csv(
        args["unmatched_peptides"], sep="\t", index=False
    )
    with open("unmatched_peptides.fasta", "w") as u:
        u.write(unmatched_peptides[1])
