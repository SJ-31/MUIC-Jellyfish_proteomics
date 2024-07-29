#!/usr/bin/env python


import pandas as pd


def removeFrom(target, string):
    target = set(target)
    if string in target:
        target.remove(string)
    return list(target)


def splitProteins(row):
    new_len = len(row["protein"])
    new = row.drop("protein").map(lambda x: [x] * new_len)
    new["protein"] = row["protein"]
    new_df = pd.DataFrame(new.to_dict())
    new_df = new_df[[new_df.columns[-1]] + new_df.columns[:-1].to_list()]
    return new_df


def cleanListCol(colname, frame, split):
    return (
        frame[colname]
        .apply(str.split, sep=split)
        .apply(removeFrom, string="")
        .apply(removeFrom, string="NA")
    )


def readDirectlfq(directlfq):
    df = pd.read_csv(directlfq, sep="\t").iloc[:, 1:]
    rename_mapping = dict(
        zip(df.columns[1:], [f"directlfq-{col}" for col in df.columns[1:]])
    )
    df = df.rename(columns=rename_mapping)
    dlfq_copy = df.copy()
    dlfq_copy["protein"] = cleanListCol("protein", dlfq_copy, ";").apply(list)
    dlfq_copy = dlfq_copy[dlfq_copy["protein"].map(lambda x: len(x) > 0)]
    split = [splitProteins(row[1]) for row in dlfq_copy.iterrows()]
    dlfq_copy = pd.concat(split)
    dlfq_copy = dlfq_copy[~dlfq_copy["protein"].str.contains("rev_")]
    dlfq_copy = dlfq_copy.groupby("protein").median()
    dlfq_copy["ProteinId"] = dlfq_copy.index
    return dlfq_copy


def readFlashlfq(file):
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
    protdf.rename({"Protein Groups": "ProteinId"}, axis="columns", inplace=True)
    return protdf


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dlfq")
    parser.add_argument("--dlfq_sorted")
    # parser.add_argument("--flfq_sorted")
    # parser.add_argument("-f", "--flfq")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    dlfq = readDirectlfq(args["dlfq"])
    dlfq.to_csv(args["dlfq_sorted"], sep="\t", na_rep="NA", index=False)
    # flfq = readFlashlfq(args["flfq"])
    # flfq.to_csv(args["flfq_sorted"], sep="\t", na_rep="NA", index=False)
