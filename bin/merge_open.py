#!/usr/bin/env python

import re
import pandas as pd


def clean_peptide(peptide):
    peptide = peptide.upper().replace("X", "")
    return "".join(re.findall("[A-Z]+", peptide))


def get_fasta_str(df, col):
    copy = df.copy()
    copy = copy[~copy["peptideIds"].isna()]
    copy["clean"] = copy["peptideIds"].apply(clean_peptide)
    fasta_str = "".join(
        copy.apply(
            lambda x: f'>{x["ProteinId"]}\n{x["clean"]}\n', axis="columns"
        ).to_list()
    )
    return fasta_str


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--intersected_searches")
    parser.add_argument("-s", "--open_searches")
    parser.add_argument("-d", "--database_output")
    parser.add_argument("--unknown_output")
    parser.add_argument("-k", "--unknown_output_fasta")
    parser.add_argument("--unmatched_peptides_in")
    parser.add_argument("--unmatched_peptides_out")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


def concat_col(df, col1, col2) -> pd.Series:
    """
    Concatenate columns col1 and col2 such that their entries are
    joined end-to-end in a csv list.
    """
    return df[col1].combine(df[col2], lambda x, y: f"{x};{y}")


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
    shared["peptideIds"] = concat_col(shared, "peptideIds_x", "peptideIds_y")
    shared["q.value"] = concat_col(shared, "q.value_x", "q.value_y")
    shared["posterior_error_prob"] = concat_col(
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
    open_searches = open_searches[~open_searches["ProteinId"].isin(shared["ProteinId"])]

    final = pd.concat([shared, open_searches, proteins])
    final["inferred_by"] = "initial_database"
    db_hits = final.query("ProteinId.str.contains('P')")
    unknown_hits = final.query("ProteinId.str.contains('O|T|D')")
    fasta_str = get_fasta_str(unknown_hits, "seq")
    return {
        "all": final,
        "database_hits": db_hits,
        "unknown_hits": unknown_hits,
        "fasta_str": fasta_str,
    }


def check_unmatched(final_frame, peps):
    """
    Filter out the peptides in the dataframe "peps" to leave only peptides that
    are truly unmatched
    """
    pep_df = pd.read_csv(peps, sep="\t")
    all_peps = pd.Series(
        [
            pep
            for pep_list in final_frame["peptideIds"].str.split(";")
            for pep in pep_list
        ]
    )
    actual_unmatched = pep_df.query("~(`peptideIds`.isin(@all_peps))")
    fasta_str = get_fasta_str(actual_unmatched, "peptideIds")
    return (actual_unmatched, fasta_str)


if __name__ == "__main__":
    args = parse_args()
    m = main(args)
    m["database_hits"].to_csv(
        args["database_output"], sep="\t", index=False, na_rep="NA"
    )
    m["unknown_hits"].to_csv(args["unknown_output"], sep="\t", index=False, na_rep="NA")
    with open(args["unknown_output_fasta"], "w") as f:
        f.write(m["fasta_str"])
    unmatched_peptides = check_unmatched(m["all"], args["unmatched_peptides_in"])
    unmatched_peptides[0].to_csv(args["unmatched_peptides_out"], sep="\t", index=False)
    with open("unmatched_peptides.fasta", "w") as u:
        u.write(unmatched_peptides[1])
