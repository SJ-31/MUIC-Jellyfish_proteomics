# Scratch session 2024-05-04-Sat
import pandas as pd
import re
import py4cytoscape as py4

wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
data = pd.read_csv(f"{wd}/results/C_indra_A/1-First_pass/C_indra_all.tsv", sep="\t")
clusters = pd.read_csv(
    f"{wd}/tests/testthat/output/cluster_go_semantic/cluster_members.tsv",
    sep="\t",
)


def flattenedJoined(
    collection: pd.Series | list | tuple, splt: str = None
) -> list[str]:
    """
    Flatten a collection of strings completely, optionally splitting apart
    entries by according to the regex `splt`
    If a series, any NA entries are ignored
    """
    if isinstance(collection, pd.Series):
        collection = collection[collection.notna()]
    if splt:
        return [x for y in collection for x in re.split(splt, y)]


go_terms = flattenedJoined(data["GO_IDs"], ";")
chosen = clusters[["id", "hc_height_0.54"]]
