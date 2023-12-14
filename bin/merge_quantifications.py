#!/usr/bin/env python

import numpy as np
import pandas as pd
# from emPAI import calculate_emPAI # Calculate once de novo peptides have been de


def clean(target, string):
    target = set(target)
    if string in target:
        target.remove(string)
    return list(target)


def read_directlfq(directlfq):
    df = pd.read_csv(directlfq, sep="\t").iloc[:, 1:]
    rename_mapping = dict(
        zip(df.columns[1:], [f"directlfq-{col}" for col in df.columns[1:]]))
    df = df.rename(columns=rename_mapping)
    return (df)


def clean_list_col(colname, frame, split):
    return (frame[colname].apply(str.split, sep=split).apply(
        clean, string="").apply(clean, string="NA"))


def sort_val(vals):
    if len(vals) > 1:
        # Sort the VALs and return the second smallest val
        # If BOTH engines have small vals for the protein,
        #   then it will be accepted
        sorted_vals = [float(vals) for vals in vals]
        sorted_vals.sort()
        return (sorted_vals[1])
    return (float(vals[0]))


def split_proteins(row):
    new_len = len(row["protein"])
    new = row.drop("protein").map(lambda x: [x] * new_len)
    new["protein"] = row["protein"]
    new_df = pd.DataFrame(new.to_dict())
    new_df = new_df[[new_df.columns[-1]] + new_df.columns[:-1].to_list()]
    return new_df


def merge_dlfq(direct_lfq, prot_df, pep_thresh, fdr):
    direct_lfq["protein"] = clean_list_col("protein", direct_lfq,
                                           ";").apply(list)
    direct_lfq = direct_lfq[direct_lfq["protein"].map(lambda x: len(x) > 0)]
    split = [split_proteins(row[1]) for row in direct_lfq.iterrows()]
    direct_lfq = pd.concat(split)
    direct_lfq = direct_lfq[~direct_lfq["protein"].str.contains("rev_")]
    direct_lfq = direct_lfq.groupby("protein").aggregate(np.median)
    merged = pd.merge(prot_df,
                      direct_lfq,
                      how="left",
                      left_on="ProteinId",
                      right_on="protein")
    merged["posterior_error_prob"] = (clean_list_col("posterior_error_prob",
                                                     merged,
                                                     ",").apply(sort_val))
    merged["q.value"] = (clean_list_col("q.value", merged,
                                        ",").apply(sort_val))
    merged = merged[merged["posterior_error_prob"] < pep_thresh]
    merged = merged[merged["q.value"] < fdr]
    return (merged)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--direct_lfq")
    parser.add_argument("-i", "--intersected_searches")
    parser.add_argument("-p", "--pep_threshold")
    parser.add_argument("-q", "--fdr_threshold")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == '__main__':
    args = parse_args()
    dlfq = read_directlfq(args["direct_lfq"])
    proteins = pd.read_csv(args["intersected_searches"], sep="\t")
    merged = merge_dlfq(dlfq, proteins, float(args["pep_threshold"]),
                        float(args["fdr_threshold"]))
    # qdf = calculate_emPAI(merged, mapping, (350, 1600))
    merged["Annotation_method"] = "initial_database"
    merged.to_csv(args["output"], sep="\t", index=False)
