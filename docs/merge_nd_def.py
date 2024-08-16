#!/usr/bin/env python
import shutil
import polars as pl
import sys
from collections import Counter
from pathlib import Path


if "Bio_SDD" in str(Path().absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
else:
    wd = "/home/shannc/workflow"

results = f"{wd}/results"

passes = ["1-First_pass", "2-Second_pass"]
d_path = f"{results}/C_indra"
nd_path = f"{results}/ND_C_indra"
nd_merged_path = f"{results}/ND_MERGED"

CONCAT = {"query_map.tsv", "unknown_hits.tsv"}

# Set up blast query files
nd_merged_seq_path = f"{nd_path}/Databases/seq-header_mappings.tsv"
default_seq_path = f"{d_path}/Databases/seq-header_mappings.tsv"
seq_path = f"{nd_merged_path}/key_seq-header_mappings.tsv"

if not Path(seq_path).exists():
    source_seq_map = pl.read_csv(default_seq_path, separator="\t").filter(
        pl.col("header") != "DECOY"
    )
    sink_seq_map = pl.read_csv(nd_merged_seq_path, separator="\t").filter(
        pl.col("header") != "DECOY"
    )
    seq_map = source_seq_map.join(sink_seq_map, on="header", suffix="_sink", how="left")
    seq_map.write_csv(seq_path, separator="\t", null_value="NA")
else:
    seq_map = pl.read_csv(seq_path, separator="\t", null_values="NA")

source_ids2sink: pl.DataFrame = seq_map.select("id", "id_sink")

source_combined_dirs = [f"{d_path}/{p}/Combined" for p in passes]
source_combined_dirs2 = [f"{nd_path}/{p}/Combined" for p in passes]
source_blast_dirs = [f"{d_path}/{p}/Unmatched/BLAST" for p in passes]
source_blast_dirs2 = [f"{nd_path}/{p}/Unmatched/BLAST" for p in passes]

sink_combined_dirs = [f"{nd_merged_path}/{p}/Combined" for p in passes]
sink_blast_dirs = [f"{nd_merged_path}/{p}/Unmatched/BLAST" for p in passes]

# Copy all blast query files over, except for the important search results

for source, source2, sink in zip(
    source_combined_dirs, source_combined_dirs2, sink_combined_dirs
):
    files = list(filter(lambda x: ".fasta" in x.name, Path(source).iterdir()))
    for f in files:
        shutil.copy(f, sink)
    for to_concat in CONCAT:
        s1 = pl.read_csv(f"{source}/{to_concat}", separator="\t")
        s2 = pl.read_csv(f"{source2}/{to_concat}", separator="\t")
        combined = pl.concat([s1, s2])
        if to_concat == "unknown_hits.tsv":
            combined = combined.unique("ProteinId")
        combined.write_csv(f"{sink}/{to_concat}", separator="\t", null_value="NA")

blast_headers = [
    "queryId",
    "subjectId",
    "sAlignStart",
    "sAlignEnd",
    "alignLen",
    "bitscore",
    "evalue",
    "pident",
    "nident",
    "nmismatch",
    "ngaps",
]
fasts_headers = [
    "queryId",
    "subjectId",
    "pident",
    "alignLen",
    "nmismatch",
    "ngaps",
    "qAlignStart",
    "qAlignEnd",
    "sAlignStart",
    "aAlignEnd",
    "evalue",
    "bitscore",
]

for source, source2, sink in zip(
    source_blast_dirs, source_blast_dirs2, sink_blast_dirs
):
    blast_files = list(filter(lambda x: ".csv" in x.name, Path(source).iterdir()))
    blast_files_2 = list(filter(lambda x: ".csv" in x.name, Path(source2).iterdir()))
    file_tracker: dict = {"blast": [], "fasts": []}
    for f in blast_files:
        if f.stat().st_size == 0:
            shutil.copy(f, sink)
            continue
        if "blast" in f.name:
            read = pl.read_csv(f, new_columns=blast_headers)
            file_tracker["blast"].append((read, f"{sink}/{f.name}"))
        elif "fasts" in f.name:
            read = pl.read_csv(f, new_columns=fasts_headers)
            file_tracker["fasts"].append((read, f"{sink}/{f.name}"))
        swapped = (
            read.join(source_ids2sink, left_on="subjectId", right_on="id")
            .drop("subjectId")
            .rename({"id_sink": "subjectId"})
            .select(read.columns)
        )
        swapped.write_csv(f"{sink}/{f.name}", include_header=False)

    # Combine results from ND run only
    to_combine: dict = {"blast": [], "fasts": []}
    for f in blast_files_2:
        if f.stat().st_size == 0:
            continue
        if "blast" in f.name:
            read = pl.read_csv(f, new_columns=blast_headers)
            to_combine["blast"].append(read)
        elif "fasts" in f.name:
            read = pl.read_csv(f, new_columns=fasts_headers)
            to_combine["fasts"].append(read)

    for t in ["blast", "fasts"]:
        prev, filename = file_tracker[t][-1]
        concatted: pl.DataFrame = pl.concat(to_combine[t] + [prev])
        concatted.write_csv(filename, include_header=False)

for p in passes:
    current = Path(f"{nd_merged_path}/{p}")
    to_delete = list(filter(lambda x: x.is_file(), current.iterdir()))
    blast_path = current.joinpath("Unmatched/BLAST")
    cov = current.joinpath(".coverage")
    if cov.exists():
        shutil.rmtree(cov)
    to_delete.extend(list(filter(lambda x: ".csv" not in x.name, blast_path.iterdir())))
    for d in to_delete:
        d.unlink()
