#!/usr/bin/env python

import pandas as pd


def clean(target, string):
    target = set(target)
    if string in target:
        target.remove(string)
    return list(target)


def clean_list_col(colname, frame, split):
    return (
        frame[colname]
        .apply(str.split, sep=split)
        .apply(clean, string="")
        .apply(clean, string="NA")
    )


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--intersected_searches")
    parser.add_argument("-s", "--open_searches")
    parser.add_argument("-u", "--unmatched_peptides")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


def concatComma(df, col1, col2) -> pd.Series:
    """
    Concatenate columns col1 and col2 such that their entries are
    joined end-to-end in a csv list.
    """
    return df[col1].combine(df[col2], lambda x, y: f"{x},{y}")


def main(args: dict):
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

    final = pd.concat([shared, open_searches, proteins])
    final["Anno_method"] = "initial_database"
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
