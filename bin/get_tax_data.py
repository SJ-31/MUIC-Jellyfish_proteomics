#!/usr/bin/env python
import ete4 as et
import sqlite3
import re
import polars as pl


def clean_name(name: str) -> str:
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


def add_taxid_col(df: pl.DataFrame, ncbi: et.NCBITaxa):
    def id_col_helper(cur_df, organism_col: str) -> pl.DataFrame:
        taxids = ncbi.get_name_translator(cur_df[organism_col])
        return cur_df.with_columns(
            TaxID=pl.col(organism_col).map_elements(
                lambda x: taxids.get(x, [-1])[0], return_dtype=pl.Int64
            )
        )

    added = id_col_helper(df, "organism")
    has_ids = added.filter(pl.col("TaxID") != -1)
    no_ids = added.filter(pl.col("TaxID") == -1).with_columns(
        organism=pl.col("organism").map_elements(clean_name, return_dtype=pl.String)
    )
    no_ids = id_col_helper(no_ids, "organism")
    added = pl.concat([has_ids, no_ids])
    return added


def get_tax_data(df: pl.DataFrame, ranks: list[str], ncbi: et.NCBITaxa):
    """Retrieve NCBI taxon ids for organisms in `df`, and parse the
    lineage data to retrieve the Kingdom, Phylum and classes of the organisms
    """
    added = add_taxid_col(df, ncbi)
    tree = ncbi.get_topology(
        filter(lambda x: x != -1, added["TaxID"]), intermediate_nodes=True
    )

    def id_lookup(id: int):
        find = tree.search_nodes(name=str(id))
        return next(find, None)

    def get_lineage(taxid: int, translate=True) -> str:
        """Filters the lineage of the specified taxid by the given ranks, returning
        a concatenated string of the lineages"""
        found = id_lookup(taxid)
        if not found:
            return ""
        found = found.props
        rank_map = {k: v for v, k in ncbi.get_rank(found["lineage"]).items()}
        filtered_ranks: list[int] = [rank_map.get(r.lower(), -1) for r in ranks]
        if not translate:
            return ";".join(filtered_ranks)
        translated = ncbi.get_taxid_translator(filtered_ranks)
        joined = ";".join([translated.get(r, "NA") for r in filtered_ranks])
        return joined

    lineages = (
        added["TaxID"]
        .map_elements(lambda x: get_lineage(x, ranks), return_dtype=pl.String)
        .str.split_exact(";", len(ranks) - 1)
        .struct.rename_fields(ranks)
        .alias("fields")
        .to_frame()
        .unnest("fields")
    )
    return pl.concat(
        [added.select("ProteinId", "TaxID"), lineages], how="horizontal"
    ).filter(pl.col("TaxID") != -1)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-n", "--ncbi_taxdump")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    RANKS: list = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    args = parse_args()
    database_locked = True
    while database_locked:
        try:
            NCBI = et.NCBITaxa(taxdump_file=args["ncbi_taxdump"])
            database_locked = False
        except Exception as e:
            if not isinstance(e, sqlite3.OperationalError):
                raise
    df = pl.read_csv(args["input"], separator="\t", null_values="NA").filter(
        pl.col("organism").is_not_null()
    )
    tax = get_tax_data(df, ranks=RANKS, ncbi=NCBI)
    tax.write_csv(args["output"], separator="\t")
