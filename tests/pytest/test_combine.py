#!/usr/bin/env ipython

import pandas as pd


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def strInSeries(series: pd.Series, query: str):
    return series.apply(lambda x: isinstance(x, float) or query in x).any()


first_pass = "./results/jellyfish/1-First_pass/"
second_pass = "./results/jellyfish/2-First_pass/"
blast_vars = ["nO_nD", "O_D", "O_nD", "nO_D"]
to_check = ["Unmatched/eggNOG", "Unmatched/InterPro", "Unmatched/BLAST"]
prefix_string = ["eggnog_meta-", "interpro_meta-", "db_hits-"]
caught = []
for b in blast_vars:
    for c, p in zip(to_check, prefix_string):
        current = f"{first_pass}/{c}/{p}{b}.tsv"
        df = read_tsv(current)
        print("___")
        print(current)
        print(df.isna().sum())
        # print(df["ID_method"].value_counts())
        if strInSeries(df["ID_method"], "-"):
            caught.append(current)
        if c == "Unmatched/BLAST":
            current2 = f"{first_pass}/{c}/{b}_unmatched.tsv"
            df2 = read_tsv(current2)
            if strInSeries(df2["ID_method"], "-"):
                caught.append(current2)
print(caught)
# idf = read_tsv(
#     f"{first_pass}/{to_check[1]}/{prefix_string[1]}{blast_vars[1]}.tsv"
# )
