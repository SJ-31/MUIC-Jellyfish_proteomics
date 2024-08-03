import polars as pl
import sys
import pathlib
import re

wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
ref = f"{wd}/data/reference"
toxins = f"{ref}/toxin_groups.tsv"
protein_groups = f"{wd}/config/protein_groups.toml"
cog_tsvs = f"{ref}/cog_ipr_pf_pthr"
cog_map = f"{ref}/cog_go_mapping.tsv"
nx_graph = f"{ref}/go_networkx.gml"
go_data = f"{ref}/with_levels.tsv"
sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
import helpers as hh
import go_subset as gs
import grouping as gg


def accession2db(acc) -> str:
    if re.match("^PTHR.*", acc):
        return "PANTHER"
    elif re.match("^IPR.*", acc):
        return "INTERPRO"
    elif re.match("^PF.*", acc):
        return "PFAM"


def db_formatting(file, group) -> pl.DataFrame:
    return (
        pl.read_csv(file, separator="\t")
        .rename({"accession": "Accession", "name": "Name"})
        .with_columns(
            Source=pl.col("Accession").map_elements(
                accession2db, return_dtype=pl.String
            ),
            Group=pl.lit(group),
        )
        .drop("type", "interpro_accession")
    )


toxin_grouping = pl.read_csv(toxins, separator="\t").with_columns(
    Group=pl.lit("venom_component")
)
db_files = pathlib.Path(cog_tsvs)
all_tsv_grouping: pl.DataFrame = pl.concat(
    [db_formatting(f, f.name.replace(".tsv", "")) for f in db_files.iterdir()]
).filter(~pl.col("Accession").is_in(toxin_grouping["Accession"]))
all_tsv_grouping = pl.concat([all_tsv_grouping, toxin_grouping])
all_tsv_grouping.write_csv(cog_map, separator="\t")
gs.go_group_mapping(nx_graph, go_data, protein_groups, cog_map, "cog")
