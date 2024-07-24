#!/usr/bin/env python

import re
import numpy as np
import polars as pl


def clean_peptide(peptide):
    if re.search("[a-z]", peptide):
        mod_regex = re.compile(r"\[[A-Za-z_]+\:[_A-Za-z]+\]")
        peptide = re.subn(mod_regex, "", peptide)[0]
    peptide = peptide.replace("X", "")
    peptide = re.sub("^n", "", peptide)
    return "".join(re.findall("[A-Z]+", peptide))


def get_fasta_str(df):
    exploded = (
        df.with_columns(pl.col("peptideIds").str.split(";"))
        .explode("peptideIds")
        .with_columns(
            pl.col("peptideIds")
            .map_elements(clean_peptide, return_dtype=pl.String)
            .alias("cleaned")
        )
    )
    other_cols = list(exploded.columns)
    other_cols.remove("ProteinId")
    query_map = (
        exploded.group_by("ProteinId")
        .agg(pl.int_range(pl.len()).alias("id_index"), pl.col(other_cols))
        .explode(["id_index"] + other_cols)
        .with_columns(queryId=pl.concat_str(["ProteinId", "id_index"], separator="."))
    ).select("queryId", "ProteinId", "cleaned")
    fastas = []
    for i in np.array_split(query_map, 5):
        df_i = pl.DataFrame(i)
        df_i.columns = query_map.columns
        fasta_str = []
        for q, seq in zip(df_i["queryId"], df_i["cleaned"]):
            fasta_str.append(f">{q}\n{seq}\n")
        fastas.append("".join(fasta_str))
    return query_map, fastas


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--intersected_searches")
    parser.add_argument("-s", "--open_searches")
    parser.add_argument("-d", "--database_output")
    parser.add_argument("-q", "--query_map_output")
    parser.add_argument("--unknown_output")
    parser.add_argument("-k", "--unknown_output_fasta")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


def main(args: dict):
    proteins = pl.read_csv(
        args["intersected_searches"], separator="\t", null_values="NA"
    )
    open_searches = pl.read_csv(
        args["open_searches"], separator="\t", null_values="NA"
    ).rename({"q-value": "q.value"})

    to_concat = [
        "peptideIds",
        "q.value",
        "posterior_error_prob",
        "ID_method",
        "ProteinGroupId",
    ]
    concat_exprs = [
        pl.concat_str([f"{i}", f"{i}_right"], separator=";").alias(i) for i in to_concat
    ]
    shared = (
        proteins.join(
            open_searches.select(["ProteinId"] + to_concat),
            on="ProteinId",
        )
        .with_columns(concat_exprs)
        .select(proteins.columns)
    )
    proteins = proteins.filter(~pl.col("ProteinId").is_in(shared["ProteinId"]))
    open_searches = open_searches.filter(
        ~pl.col("ProteinId").is_in(shared["ProteinId"])
    )
    final = pl.concat([proteins, shared, open_searches]).with_columns(
        inferred_by=pl.lit("initial_database")
    )
    db_hits = final.filter(pl.col("ProteinId").str.contains("P"))
    unknown_hits = final.filter(pl.col("ProteinId").str.contains("T|D"))
    query_map, fasta_strs = get_fasta_str(unknown_hits)
    return {
        "all": final,
        "database_hits": db_hits,
        "unknown_hits": unknown_hits,
        "fasta_strs": fasta_strs,
        "query_map": query_map,
    }


if __name__ == "__main__":
    args = parse_args()
    m = main(args)
    m["database_hits"].write_csv(
        args["database_output"], separator="\t", null_value="NA"
    )
    m["unknown_hits"].write_csv(args["unknown_output"], separator="\t", null_value="NA")
    m["query_map"].write_csv(args["query_map_output"], separator="\t", null_value="NA")
    for i, s in enumerate(m["fasta_strs"]):
        name = args["unknown_output_fasta"].replace(".fasta", f"_{i}.fasta")
        with open(name, "w") as f:
            f.write(s)
