#!/usr/bin/env python
#
import re
import pandas as pd
import numpy as np


def write_unmatched(blast_df, failed_filter, prot_df, tsv_name):
    # Write the entries unmatched by blast to a new fasta and tsv file
    # all_df = pd.concat([prot_df, failed_filter])
    unmatched_fasta = ""
    unmatched = (prot_df[~prot_df["ProteinId"].isin(blast_df["queryID"])])
    unmatched = pd.concat([prot_df, failed_filter])
    unmatched.to_csv(tsv_name, sep="\t", index=False)
    for row in unmatched.iterrows():
        entry = f">{row[1]['ProteinId']}\n{row[1]['seq']}\n"
        unmatched_fasta = unmatched_fasta + entry
    return unmatched_fasta


def clean_peptide(peptide):
    peptide = peptide.upper().replace("X", "")
    return ''.join(re.findall("[A-Z]+", peptide))


def group_peps(df):
    '''
    Concatenate values in rows, separating by commas
    '''
    grouped = {}
    df["peptideIds"] = [clean_peptide(pep) for pep in df["peptideIds"]]
    for col in df.columns:
        joined_up = set([c for c in df[col]])
        if "lfq" in col:
            joined_up = [np.mean(list(joined_up))]
        grouped[col] = ",".join([str(j) for j in joined_up])
    return pd.DataFrame(grouped, index=[0])


def adjust_pep(df):
    n_matched = df.shape[0]
    df["posterior_error_prob"] = df["posterior_error_prob"].astype(
        float) * n_matched
    return df


def merge_blast(b_df, prot_df, keep_best_only, ident_thresh, e_thresh,
                pep_thresh, one_hit_wonders):
    '''
    Optionally, filter blast results by percent identity and evalue, then
    merge with the protein identifications dataframe
    '''
    b_df = b_df[(b_df["pident"] >= ident_thresh)
                & (b_df["evalue"] <= e_thresh)]
    b_df = b_df.sort_values(by="evalue", kind="stable")
    joined = pd.merge(b_df, prot_df, left_on="queryID", right_on="ProteinId")
    copy = joined.copy()
    if keep_best_only:
        # Keep only the blast hit with the lowest e-value
        # If this is False, then degenerate peptides are allowed
        joined = joined.groupby("queryID").nth(1).reset_index()
    else:
        # Adjust Percolator PEPs by the number of proteins at peptide
        # is matched with
        joined = joined.groupby("ProteinId").apply(adjust_pep)
        joined = joined[joined["posterior_error_prob"] <= pep_thresh]
    did_not_pass = (copy[~copy["queryID"].isin(joined["queryID"])].groupby(
        "queryID").nth(1)[["ProteinId", "seq"]])
    # Extract the psms that failed any filters
    joined = (joined.groupby(["subjectID"
                              ]).apply(group_peps).reset_index(drop=True))
    if not one_hit_wonders:
        # Filter out protein identifications that are matched by only
        # one peptide
        return (joined[joined["queryID"].str.contains(",")], did_not_pass)
    return (joined, did_not_pass)


def known_from_database(blast_df, db_df):
    '''
    Append de novo peptides already found in database search properly back into
    the search output
    '''
    already_found = pd.merge(db_df,
                             blast_df.filter(["subjectID", "seq"]),
                             left_on="ProteinId",
                             right_on="subjectID",
                             how="left")
    already_found["peptideIds"] = (already_found["peptideIds"].str.cat(
        already_found["seq_y"].to_list(), sep=","))
    already_found = (already_found.drop(
        ["seq_y", "subjectID"], axis="columns").rename({"seq_x": "seq"},
                                                       axis="columns"))
    return already_found


def blast_id_only(blast_df, db_df, mapping_df):
    '''
    Find proteins identified solely from blast hits and make
    a identification dataframe out of them
    '''
    not_found = blast_df[~blast_df["subjectID"].isin(db_df["ProteinId"])]
    blast_only = pd.merge(not_found,
                          mapping_df,
                          left_on="subjectID",
                          right_on="id")
    blast_only["ProteinId"] = blast_only["id"]
    blast_only["seq"] = blast_only["seq_y"]
    blast_only["mass"] = blast_only["mass_y"]
    blast_only["length"] = blast_only["length_y"]
    blast_only["header"] = blast_only["header_y"]
    blast_only = (
        blast_only.loc[:,
                       ~blast_only.columns.str.contains("[xy]", regex=True)].
        loc[:, db_df.columns])
    return blast_only


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--blast_path")
    parser.add_argument("-u", "--unknown_hits")
    parser.add_argument("-d", "--database_hits")
    parser.add_argument("-m", "--mapping")
    parser.add_argument("-f", "--unmatched_fasta")
    parser.add_argument("-t", "--unmatched_tsv")
    parser.add_argument("-o", "--output")
    parser.add_argument("-i", "--identity_threshold")
    parser.add_argument("-p", "--pep_threshold")
    parser.add_argument("-e", "--evalue_threshold")
    parser.add_argument("--keep_best")
    parser.add_argument("--one_hit")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == '__main__':
    args = parse_args()
    blast_df = pd.read_csv(args["blast_path"]).dropna(axis="columns")
    mapping = pd.read_csv(args["mapping"], sep="\t")
    unknown_df = pd.read_csv(args["unknown_hits"], sep="\t")
    group_df = pd.read_csv(args["database_hits"], sep="\t")
    joined = merge_blast(blast_df,
                         unknown_df,
                         keep_best_only=bool(int(args["keep_best"])),
                         one_hit_wonders=bool(int(args["one_hit"])),
                         ident_thresh=float(args["identity_threshold"]),
                         pep_thresh=float(args["pep_threshold"]),
                         e_thresh=float(args["evalue_threshold"]))
    unmatched = write_unmatched(blast_df, joined[1], unknown_df,
                                args["unmatched_tsv"])
    with open(args["unmatched_fasta"], "w") as f:
        f.write(unmatched)
    in_db = known_from_database(joined[0], group_df)
    from_blast = blast_id_only(joined[0], group_df, mapping)
    final = pd.concat([in_db, from_blast])
    final["Annotation_method"] = "blast"
    final.to_csv(args["output"], sep="\t", index=False)

# test command
# %run ./sort_blast.py -b ../results/jellyfish/1-First_pass/unknown-blast.csv -u ../results/jellyfish/1-First_pass/Combined/unknown_hits.tsv -m ../results/jellyfish/Databases/seq-header_mappings.tsv -f "../tests/blast/unmatched.fasta" -o "../tests/blast/matched.tsv" -d ../results/jellyfish/1-First_pass/Combined/database_hits.tsv -i 70 -e 100 --keep_best 0 --one_hit 0 -p 0.05
