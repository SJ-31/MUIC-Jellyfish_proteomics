#!/usr/bin/env python

from pathlib import Path
import pandas as pd

prefix = "C_indra"
runs = ("1-First_pass", "2-Second_pass")
basedir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
found = {}
output = {
    "interpro": f"{prefix}_interpro_matched",
    "eggnog": f"{prefix}_eggnog_matched",
    "downloads": f"{prefix}_downloads_anno-3",
    "combined": f"{prefix}_all_wcoverage",
}

# Get result files, putting them in "found"
for run in runs:
    path = f"{basedir}/results/{prefix}/{run}/"
    results_path = Path(path)
    for name, file in output.items():
        try:
            found[f"{run}_{name}"] = pd.read_csv(
                list(results_path.rglob(f"{file}.tsv"))[0], sep="\t"
            )
        except IndexError:
            print(f"Not found: {run}/{file}")


def mySearch(regex: str, df, case=False):
    """Search all the text columns of `df`, return rows with any matches."""
    textlikes = df.select_dtypes(include=[object, "string"])
    return df[
        textlikes.apply(
            lambda column: column.str.contains(
                regex, regex=True, case=case, na=False
            )
        ).any(axis=1)
    ]


def myColumnsWith(regex: str, df, case=False):
    """Return columns of `df` that contain the regex"""
    contained = []
    for col in df.select_dtypes(include=[object, "string"]):
        try:
            if (
                df[col]
                    .apply(str)
                    .str.contains(regex, regex=True, case=case, na=False)
                    .any()
            ):
                contained.append(col)
        except AttributeError:
            from IPython.core.debugger import set_trace

            set_trace()  # fmt: skip
    return contained


def hasDuplicates(df, col: str) -> bool:
    if len(df[col]) != len(df[col].unique()):
        return True
    return False


def findDuplicates(df, col: str) -> pd.Series:
    return df[col][df[col].duplicated()]


duplicates_proteinid = {}
col = "ProteinId"
for name, f in found.items():
    if hasDuplicates(f, col):
        duplicates_proteinid[name] = f[f[col].isin(findDuplicates(f, col))]
        print(f"Has duplicates: {name}")
