#!/usr/bin/env python

import polars.selectors as cs
import polars as pl
import pathlib


ANNO_COLS = [
    "PANTHER",
    "seed_ortholog",
    "COG_category",
    "KEGG_Genes",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "interpro_accession",
    "interpro_description",
    "interpro_pathways",
    "interpro_db",
    "GO",
    "GO_evidence",
    "PFAMs",
    "eggNOG_OGs",
    "Description",
    "Preferred_name",
]
ID_COLS = [
    "header",
    "organism",
    "lineage",
    "NCBI_ID",
    "UniProtKB_ID",
    "inferred_by",
]


def save_seen_anno(
    rpath: str, prefix: str, output: str, write: bool = False
) -> pl.DataFrame:
    try:
        previous_saved = pl.read_csv(output, separator="\t", null_values="NA")
    except FileNotFoundError:
        previous_saved = None
    to_read = []
    for p in ["1-First_pass", "2-Second_pass"]:
        to_read.extend(
            [
                f"{rpath}/{p}/Unmatched/Database-annotated/{prefix}_downloads_anno-complete.tsv",
                f"{rpath}/{p}/Unmatched/eggNOG/{prefix}_eggnog_matched.tsv",
                f"{rpath}/{p}/Unmatched/InterPro/{prefix}_interpro_matched.tsv",
            ]
        )
    all_annotated = [
        pl.read_csv(a, separator="\t", null_values="NA", ignore_errors=True)
        for a in to_read
    ]
    together = (
        pl.concat(all_annotated, how="diagonal_relaxed")
        .select(ID_COLS + ANNO_COLS)
        .with_columns(
            pl.col("*").replace_strict(old="-", new=None, default=pl.col("*"))
        )
    )
    if previous_saved:
        together = together.filter(~pl.col("header").is_in(previous_saved["header"]))
        together = pl.concat([together, previous_saved])
    mask = together.with_columns(pl.col(ANNO_COLS).is_not_null()).with_columns(
        should_keep=pl.any_horizontal(ANNO_COLS)
    )
    together = together.filter(mask["should_keep"])
    uniques = together.unique("header")
    if write:
        uniques.write_csv(output, separator="\t", null_value="NA")
    return uniques


def save_seen_eggnog(rpath: str, header_mapping: str, write_to: str = ""):
    header_map = pl.read_csv(header_mapping, separator="\t")
    to_read = [
        f"{rpath}/{p}/Unmatched/eggNOG" for p in ["1-First_pass", "2-Second_pass"]
    ]

    def helper(path):
        eg = pathlib.Path(path)
        annotations = next(eg.glob("*annotations"))
        orthologs = next(eg.glob("*.orthologs"))
        seed_orthologs = next(eg.glob("*seed_orthologs"))
        result = []
        for file, qcol in zip(
            [annotations, orthologs, seed_orthologs], ["#query", "#query", "#qseqid"]
        ):
            df = (
                pl.read_csv(file, separator="\t", comment_prefix="##", null_values="NA")
                .join(header_map, left_on=qcol, right_on="id")
                .drop([qcol, "seq", "mass", "length"])
            )
            result.append(df)
        return result

    first = helper(to_read[0])
    second = helper(to_read[1])
    results = []
    for i in range(3):
        df = pl.concat([first[i], second[i]]).unique("header")
        results.append(df)

    if write_to:
        for df, output in zip(results, ["annotations", "orthologs", "seed_orthologs"]):
            file = f"{write_to}/eggnog_{output}.tsv"
            df.write_csv(file, separator="\t", null_value="NA")
    return results


c_indra_a = "./results/C_indra_A"

prefix = "C_indra"
results = f"./results/{prefix}"
header_mapping = f"{results}/Databases/seq-header_mappings.tsv"

# Output directory to put the saved files into
saved = f"{results}/.saved"
annotations_path = f"{saved}/annotations.tsv"
eggnog_path = f"{saved}/eggnog.tsv"
interpro_path = f"{saved}/interpro.tsv"

# anno: pl.DataFrame = save_seen_anno(c_indra_a, "C_indra", annotations, write=True)
# eggnog = save_seen_eggnog(results, header_mapping, "./results/C_indra_A/.saved")

hh.retrieve_saved(
    "./query.fasta",
    "./results/C_indra_A/.saved",
    "./results/C_indra_A/Databases/seq-header_mappings.tsv",
    "eggnog",
)
