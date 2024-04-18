import sys

import pandas as pd
from kegg_pull import map as kpm
import re

from pandas._libs.missing import NAType

ALL_MAPPINGS: dict = {}  # Dictionary containing dictonaries of all
# organism-specific mappings


def flattenRec(collection) -> list:
    """
    Recursively flatten nested lists
    """
    if not collection:
        return []
    current = collection[0]
    if isinstance(current, set):
        current = list(current)
    if isinstance(current, list):
        return flattenRec(current) + flattenRec(collection[1:])
    return collection[:1] + flattenRec(collection[1:])


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


def gene2pathway(
    kegg_gene: str, all_mappings: dict, sep: str = ";"
) -> NAType | str:
    org: str = re.sub(":.*", "", kegg_gene)
    mapping: dict = all_mappings.get(org)
    if not mapping:
        mapping = kpm.database_link(org, "pathway")
        all_mappings[org] = mapping
    mapped = mapping.get(kegg_gene)
    if mapped:
        mapped = {re.sub("path:", "", m) for m in mapped}
        return sep.join(mapped)
    return pd.NA


def removeDuplicates(string) -> NAType | str:
    if pd.isna(string):
        return string
    splits = set(re.split("[;,]", string))
    return ";".join(splits)


def keggItemFromGenes(gene_str, db) -> NAType | str:
    items = set()
    splits = re.split("[;,]", gene_str)
    match db:
        case "ko":
            get_rid = "ko:"
        case "enzyme":
            get_rid = "ec:"
        case "module":
            get_rid = "md:"
        case _:
            get_rid = ""
    if db == "pathway":
        for gene in splits:
            items.add(gene2pathway(gene, ALL_MAPPINGS))
        if pd.NA in items:
            items.remove(pd.NA)
    else:
        for s in splits:
            try:
                found: dict = kpm.entries_link([s], db)
                for query, result in found.items():
                    for r in result:
                        items.add(re.sub(get_rid, "", r))
            except RuntimeError:
                print(f"WARNING: KEGG item {s} failed mapping")
                continue
    if not items:
        return pd.NA
    return ";".join(items)


def joinNotNa(x, y) -> str | NAType:
    if not pd.isna(x) and not pd.isna(y):
        return f"{x};{y}"
    elif pd.isna(x) and pd.isna(y):
        return pd.NA
    elif pd.isna(x):
        return y
    return x


def mapInDf(df, target_db, target_col) -> pd.DataFrame:
    """
    :param df:
    :param target_db:  One of "pathway" or "module"
    :param target_col:  One of "KEGG_Pathway" or "KEGG_Module"
    :return:
    """
    df["KEGG_Genes"] = df["KEGG_Genes"].apply(removeDuplicates)
    has_kegg = df[~df["KEGG_Genes"].isna()]
    mapped_pathways = has_kegg["KEGG_Genes"].apply(
        lambda x: keggItemFromGenes(x, target_db)
    )
    has_kegg[target_col] = has_kegg[target_col].combine(
        mapped_pathways, lambda x, y: joinNotNa(x, y)
    )
    no_kegg = df[df["KEGG_Genes"].isna()]
    joined = pd.concat([has_kegg, no_kegg])
    joined[target_col] = joined[target_col].apply(removeDuplicates)
    return joined


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__" and len(sys.argv) > 1:
    args = parse_args()
    df = pd.read_csv(args["input"], sep="\t")
    mapped = mapInDf(df, "pathway", "KEGG_Pathway")
    mapped.to_csv(args["output"], sep="\t", index=False)

