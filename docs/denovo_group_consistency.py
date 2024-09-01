import pandas as pd
import re
import sys
import polars as pl
from pathlib import Path

if "Bio_SDD" in str(Path().absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
else:
    wd = "/home/shannc/workflow"


outdir = f"{wd}/docs/figures/denovo_matches"
python_source = f"{wd}/bin"
sys.path.append(python_source)
import helpers as hh

passes = ["1-First_pass", "2-Second_pass"]
prefixes = {
    "ND": "ND_C_indra",
    "default": "C_indra",
    "Calibrated": "C_indra.calibrated",
}
results = f"{wd}/results"

path = f"{results}/{prefixes['Calibrated']}/{passes[1]}"

results_file = f"{path}/{prefixes['Calibrated']}_all_wcoverage.tsv"
hits_file = f"{path}/Combined/intersected_searches.tsv"


df = pl.read_csv(results_file, separator="\t", null_values="NA")
hits = pl.read_csv(hits_file, separator="\t", null_values="NA")

PID: str = "peptideIds"
mp: str = "MatchedPeptideIds"
exploded = (
    df.filter(pl.col(mp).is_not_null())
    .select("ProteinId", mp)
    .with_columns(pl.col(mp).str.split(";"))
    .explode(mp)
)
denovo2prot = exploded.select("ProteinId", mp)


def group_by_unique_peptides(df: pl.DataFrame) -> pl.DataFrame:
    cleaned: pl.DataFrame = df.filter(pl.col(PID).is_not_null()).with_columns(
        pl.col(PID)
        .str.split(";")
        .map_elements(
            lambda x: list(set(map(hh.clean_peptide, x))),
            return_dtype=pl.List(pl.String),
        )
        .alias("cleaned")
    )
    return cleaned.group_by("cleaned").agg(ids=pl.col("ProteinId"))


def get_id2group(df: pl.DataFrame) -> dict:
    g = group_by_unique_peptides(df).with_row_index(name="GroupUP").explode("ids")
    return dict(zip(g["ids"], g["GroupUP"]))


id2group = get_id2group(hits)


id2group_real = get_id2group(df)

agg_real = df.group_by("GroupUP").agg(pl.col("ProteinId"))

wrong_group = []

# %%
wrong = 0
for group in agg_real["ProteinId"]:
    cur_group = 0
    for index, g in enumerate(group):
        if index == 0:
            cur_group = id2group_real[g]
        else:
            if id2group_real[g] != cur_group:
                wrong += 1

total_hits = 0
for denovo, prot in denovo2prot.iter_rows():
    if id2group.get(denovo) == id2group.get(prot):
        total_hits += 1

total_hits / len(denovo2prot)


def get_longest(seqs: list[str]) -> str:
    return sorted(seqs, key=lambda x: len(x), reverse=True)[0]


grouped = df.group_by("GroupUP").agg(
    pl.col("header").first(), pl.col("entry_name").str.join(";"), pl.col("seq")
)


ompa_g = grouped.filter(
    pl.col("entry_name")
    .str.to_lowercase()
    .str.contains_any(["outer membrane protein a", "ompa family"])
).with_columns(pl.col("seq").map_elements(get_longest, return_dtype=pl.String))


hh.write_fasta(
    ompa_g["header"],
    ompa_g["seq"],
    "/home/shannc/Downloads/thesis_testzone/ompa_g.fasta",
)


ompa = df.filter(pl.col("entry_name").str.to_lowercase().str.contains("outer membrane"))

hh.write_fasta(
    ompa["header"], ompa["seq"], "/home/shannc/Downloads/thesis_testzone/ompa.fasta"
)
