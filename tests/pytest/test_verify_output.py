#!/usr/bin/env python

from pathlib import Path
import pandas as pd

prefix = "jellyfish"
path = f"./results/{prefix}/1-First_pass/"
found = {}
output = {
    "interpro": f"{prefix}_interpro_matched",
    "eggnog": f"{prefix}_eggnog_matched",
    "downloads": f"{prefix}_downloads_anno-3",
    "combined": f"{prefix}_all",
}
results_path = Path(path)
for name, file in output.items():
    found[name] = pd.read_csv(
        list(results_path.rglob(f"{file}.tsv"))[0], sep="\t"
    )


def search(regex: str, df, case=False):
    """Search all the text columns of `df`, return rows with any matches."""
    textlikes = df.select_dtypes(include=[object, "string"])
    return df[
        textlikes.apply(
            lambda column: column.str.contains(
                regex, regex=True, case=case, na=False
            )
        ).any(axis=1)
    ]


def columnsWith(regex: str, df, case=False):
    """Return columns of `df` that contain the regex"""
    contained = []
    for col in df.select_dtypes(include=[object, "string"]):
        if df[col].str.contains(regex, regex=True, case=case, na=False).any():
            contained.append(col)
    return contained
