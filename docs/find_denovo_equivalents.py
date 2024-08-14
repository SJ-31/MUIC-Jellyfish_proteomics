import polars as pl
import pandas as pd
import polars.selectors as cs
import sys
from rapidfuzz import distance as di
from pathlib import Path

if "Bio_SDD" in str(Path().absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
else:
    wd = "/home/shannc/workflow"


outdir = f"{wd}/docs/figures/denovo_matches"
python_source = f"{wd}/bin"
sys.path.append(python_source)
import helpers as hh


chosen_pass = "2-Second_pass"
prefixes = {"ND": "ND_C_indra", "default": "C_indra"}
results = f"{wd}/results"

dfs: dict = {}
peptides: dict = {}
hits: dict = {}
perc: dict = {}
perc2: dict = {}

for p, v in prefixes.items():
    path = f"{results}/{v}/{chosen_pass}"
    dfs[p] = pl.read_csv(
        f"{path}/{v}_all_wcoverage.tsv",
        separator="\t",
        null_values="NA",
    )
    peptides[p] = set(hh.flatten_by(dfs[p]["unique_peptides"]))

was_matched: set = set(
    hh.flatten_by(
        dfs["default"].filter(pl.col("MatchedPeptideIds").is_not_null())[
            "MatchedPeptideIds"
        ]
    )
)


perc_prot = (
    pl.read_csv(
        f"{results}/{prefixes['default']}/{chosen_pass}/percolator_all.tsv",
        separator="\t",
        null_values="NA",
    )
    .with_columns(pl.col("peptideIds").str.split(";"))
    .explode("peptideIds")
).with_columns(
    pl.col("peptideIds").map_elements(hh.clean_peptide, return_dtype=pl.String)
)
perc_peps = pl.read_csv(
    f"{results}/{prefixes['default']}/{chosen_pass}/percolator_peptide_map.tsv",
    separator="\t",
    null_values="NA",
).select(perc_prot.columns)

perc: pl.DataFrame = pl.concat([perc_peps, perc_prot])
perc.write_csv(f"{outdir}/filtered_percolator.tsv")

original_denovo: pl.DataFrame = hh.fasta2df(
    f"{wd}/data/protein_databases/denovo_normal_mgf/denovo_all.fasta"
).with_columns(pl.col("seq").map_elements(str, return_dtype=pl.String))

denovo: pl.DataFrame = (perc_prot.filter(pl.col("ProteinId").str.contains("D"))).unique(
    "peptideIds"
)

default_engine_peps = set(
    hh.flatten_by(
        dfs["default"].filter(pl.col("MatchedPeptideIds").is_null())["unique_peptides"]
    )
)

number_matched_engines = len(peptides["ND"] & default_engine_peps)
peptides["ND"] = (
    peptides["ND"] - default_engine_peps
)  # Get rid of all the normal engine matches, so only potential de novo matches remain


types = ["MATCHED", "ALL"]
# Difference is that "MATCHED" uses peptides that were actually identified by engines

seq_header_map = (
    pl.read_csv(
        f"{results}/C_indra/Databases/seq-header_mappings.tsv",
        separator="\t",
        null_values="NA",
    )
    .filter(pl.col("header").str.contains("DENOVO"))
    .select("id", "seq", "header")
)

denovo_peps: list[set] = [set(denovo["peptideIds"]), set(original_denovo["seq"])]

all_denovo = pl.concat(
    [
        (
            denovo.rename({"peptideIds": "seq"})
            .join(
                seq_header_map.drop("seq"),
                left_on="ProteinId",
                right_on="id",
            )
            .select(original_denovo.columns + ["ProteinId"])
            .with_columns(from_engine_ids=pl.lit(True))
        ),
        original_denovo.with_columns(ProteinId=None, from_engine_ids=pl.lit(False)),
    ]
)

# all_denovo_peps: pl.DataFrame = de
for type, d in zip(types, denovo_peps):
    # Direct instances of engine peptides being matched to de novo peptides
    direct_matches = peptides["ND"] & d
    hh.py_cat(
        [
            "Number of direct matches between ND engine peptides and de novo peptides",
            len(direct_matches),
            f"Proportion: {len(direct_matches)/len(peptides['ND'])}",
            f"Number of total remaining ND engine peptides: {len(peptides['ND'])}",
            f"Number of matched ND engine peptides: {number_matched_engines}",
        ],
        f"{outdir}/direct_matches-{type}.txt",
    )

result_file = f"{outdir}/denovo_ND_hits-COMPLETE.tsv"
if not Path(result_file).exists():
    # Attempt to find all matched ND engine peptides in a de novo peptide
    find_in_denovo: pl.DataFrame = pl.from_pandas(
        hh.find_matches(peptides["ND"], set(all_denovo["seq"]))
    )
    find_in_denovo.write_csv(result_file, separator="\t", null_value="NA")
else:
    hits = pl.read_csv(result_file, separator="\t", null_values="NA")


df: pl.DataFrame = hits.join(
    all_denovo, left_on="best_hit", right_on="seq"
).with_columns(
    pl.struct(["query", "best_hit"])
    .map_elements(
        lambda x: di.Levenshtein.distance(x["query"], x["best_hit"]),
        return_dtype=pl.Int64,
    )
    .alias("distance")
)

df.write_csv(f"{outdir}/COMPLETE_final.tsv", separator="\t", null_value="NA")


transcriptome: pl.DataFrame = perc.filter(pl.col("ProteinId").str.contains("T"))
not_confident: pl.DataFrame = hits.filter(pl.col("similarity") < 0.7)
remaining_nd = set(not_confident["query"])
direct_matches_t = remaining_nd & set(transcriptome["peptideIds"])

hh.py_cat(
    [
        "Number of direct matches between ND engine peptides and transcriptome peptides",
        len(direct_matches_t),
        f"Proportion: {len(direct_matches_t)/len(remaining_nd)}",
        f"Number of total remaining ND engine peptides: {len(remaining_nd)}",
    ],
    f"{outdir}/direct_matches-TRANSCRIPTOME.txt",
)
# Attempt to check if remaining unconfident de novo peptides are from transcriptome
# peptides
result_file = f"{outdir}/denovo_ND_hits-TRANSCRIPTOME.tsv"
if not Path(result_file).exists():
    find_in_denovo_t: pl.DataFrame = pl.from_pandas(
        hh.find_matches(remaining_nd, set(transcriptome["peptideIds"]))
    )
    find_in_denovo_t.write_csv(result_file, separator="\t", null_value="NA")
else:
    hits_t = pl.read_csv(result_file, separator="\t", null_values="NA")
