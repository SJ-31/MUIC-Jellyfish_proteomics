#!/usr/bin/env python
import polars as pl
import sys
from collections import Counter
from pathlib import Path


if "Bio_SDD" in str(Path().absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
    taxdump = "/home/shannc/Bio_SDD/tools/taxdb/taxdump.tar.gz"
else:
    wd = "/home/shannc/workflow"
    taxdump = f"{wd}/data/reference/taxdump.tar.gz"


sys.path.append(f"{wd}/bin")
import tree_viz as tv


chosen_pass = "2-Second_pass"
prefixes = {
    "ND": "ND_C_indra",
    "default": "C_indra",
    "Calibrated": "C_indra.calibrated",
}
results = f"{wd}/results"
chosen_prefix = prefixes["Calibrated"]
tax_path = f"{results}/{chosen_prefix}/{chosen_pass}/{chosen_prefix}_taxonomy.tsv"
data_path = f"{results}/{chosen_prefix}/{chosen_pass}/{chosen_prefix}_all_wcoverage.tsv"
outdir = f"{wd}/docs/figures/taxonomy"
grouped_tax_path = f"{outdir}/grouped_taxa.tsv"

data = pl.read_csv(data_path, separator="\t", null_values="NA").select(
    "ProteinId", "GroupUP"
)
tax = pl.read_csv(tax_path, separator="\t", null_values="NA")
wanted_cols = tax.columns
wanted_cols.remove("ProteinId")


def get_mode(items: list | pl.Series):
    if len(items) > 0:
        return sorted(Counter(items).items(), key=lambda x: x[1], reverse=True)[0][0]


grouped_tax: pl.DataFrame = (
    tax.join(data, on="ProteinId")
    .group_by("GroupUP")
    .agg(wanted_cols)
    .with_columns(pl.col(wanted_cols).map_elements(get_mode, return_dtype=pl.String))
)
grouped_tax.write_csv(grouped_tax_path, separator="\t", null_value="NA")


if not Path(outdir).exists():
    Path(outdir).mkdir()

tree = tv.TaxaTree(grouped_tax_path, taxdump)
phyla = tree.get_subtree(rank="phylum")
cnidaria = tree.get_subtree(sci_name="Cnidaria", rank="order")
tv.show(
    phyla,
    legendFun=lambda x: tv.rank_legend(x, "phylum"),
    save_to=f"{outdir}/Phyla.png",
    save_params={"w": 1800, "h": 1200},
)
tv.show(
    cnidaria,
    legendFun=lambda x: tv.rank_legend(x, "order", ("kingdom",)),
    save_to=f"{outdir}/cnidaria.png",
    save_params={"w": 1800, "h": 1200},
)
