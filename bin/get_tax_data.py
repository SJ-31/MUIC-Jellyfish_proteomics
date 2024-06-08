#!/usr/bin/env python
import pytaxonkit as pk
import re
import pandas as pd
import polars as pl


def cleanName(name: str) -> str:
    if " sp" in name or " subsp" in name:
        name = re.sub(r"(sub)*sp\.*.*", "", name)
    if "complex" in name:
        name = re.sub(r"complex.*", "", name)
    if "serotype" in name:
        name = re.sub(r"serotype.*", "", name)
    if " bv." in name:
        name = re.sub(r" bv\..*", "", name)
    if "(strain" in name:
        name = re.sub(r"\(strain.*", "", name)
    if "uncultured" in name:
        name = name.replace("uncultured", "")
    return name.strip()


def addTaxIdCol(df: pl.DataFrame) -> pl.DataFrame:
    id_df = pk.name2taxid(df["organism"].unique())
    failed = id_df.query("TaxID.isna()")
    cleaned_names = failed["Name"].apply(cleanName)
    try_cleaned = pk.name2taxid(cleaned_names)
    id_df["TaxID"] = id_df["TaxID"].combine(
        try_cleaned["TaxID"], lambda x, y: x if pd.isna(y) else y
    )
    taxids: dict = dict(zip(id_df["Name"], id_df["TaxID"]))
    df = df.with_columns(
        TaxID=pl.col("organism").map_elements(
            lambda x: taxids.get(x), return_dtype=pl.Int64
        )
    )
    return df


RANKS: list = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]


def getTaxData(df: pd.DataFrame) -> pl.DataFrame:
    """Retrieve NCBI taxon ids for organisms in `df`, and parse the
    lineage data to retrieve the Kingdom, Phylum and classes of the organisms
    """
    fstr = "{k};{p};{c};{o};{f};{g}"
    df = addTaxIdCol(df).filter(pl.col("TaxID").is_not_null())
    lineages = pl.from_dataframe(pk.lineage(df["TaxID"], formatstr=fstr))
    kpc = (
        lineages["Lineage"]
        .str.split_exact(";", 5)
        .struct.rename_fields(RANKS)
        .alias("fields")
        .to_frame()
        .unnest("fields")
    )
    kpc = pl.concat([lineages.select("TaxID", "Name", "Rank"), kpc], how="horizontal")
    return (
        df.join(
            kpc,
            on="TaxID",
            how="left",
        )
        .select("ProteinId", "TaxID", "Rank", *RANKS)
        .unique()
    )


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parseArgs()
    df = pl.read_csv(args["input"], separator="\t", null_values="NA")
    tax = getTaxData(df)
    tax.write_csv(args["output"], separator="\t")
