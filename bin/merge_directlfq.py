import numpy as np
from re import split
import pandas as pd

directlfq = "results/benchmark_Ecoli/1-First_pass/1-Quantify/directlfq_prot.tsv"
prot_list = "results/benchmark_Ecoli/1-First_pass/1-Combined/intersected_searches.tsv"
pep_thresh = 0.05

df = pd.read_csv(directlfq, sep="\t").iloc[:, 1:]
prot_df = pd.read_csv(prot_list, sep="\t")
rename_mapping = dict(zip(df.columns[1:],
                          [f"directlfq-{col}" for col in df.columns[1:]]))
df = df.rename(rename_mapping)

def clean(target, string):
    target = set(target)
    if string in target:
        target.remove(string)
    return list(target)

def clean_list_col(colname, frame, split):
    return (frame[colname].apply(str.split, sep=split)
            .apply(clean, string="")
            .apply(clean, string="NA"))

def sort_pep(peps):
    if len(peps) > 1:
    # Sort the PEPs and return the second smallest pep
    # If BOTH engines have small peps for the protein,
    #   then it will be accepted
        peps = [float(peps) for peps in peps]
        peps.sort()
        return(peps[1])
    return(float(peps[0]))

def split_proteins(row):
    new_len = len(row["protein"])
    new = row.drop("protein").map(lambda x: [x] * new_len)
    new["protein"] = row["protein"]
    new_df = pd.DataFrame(new.to_dict())
    new_df = new_df[[new_df.columns[-1]] + new_df.columns[:-1].to_list()]
    return new_df

all_dfs = pd.DataFrame()
df["protein"] = clean_list_col("protein", df, ";").apply(list)
df = df[df["protein"].map(lambda x: len(x) > 0)]



split = [split_proteins(row[1]) for row in df.iterrows()]
df = pd.concat(split)
df = df[~df["protein"].str.contains("rev_")]

df = df.groupby("protein").aggregate(np.median)
merged = pd.merge(prot_df, df, how="left", left_on="ProteinId", right_on="protein")

merged["posterior_error_prob"] = (clean_list_col("posterior_error_prob",
                                                 merged, ",")
                                  .apply(sort_pep))

merged = merged[merged["posterior_error_prob"] < pep_thresh]
q()
