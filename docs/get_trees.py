#!/usr/bin/env python
import sys
from pathlib import Path

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
import tree_viz as tv

taxdump = "/home/shannc/Bio_SDD/tools/taxdb/taxdump.tar.gz"
tax = "../results/C_indra/1-First_pass/C_indra_taxonomy.tsv"
outdir = "../results/C_indra/Analysis/Taxonomy"

if not Path(outdir).exists():
    Path(outdir).mkdir()

tree = tv.TaxaTree(tax, taxdump)
phyla = tree.getSubtree(rank="phylum")
cnidaria = tree.getSubtree(sci_name="Cnidaria", rank="order")
tv.show(
    phyla,
    legendFun=lambda x: tv.rankLegend(x, "phylum"),
    save_to=f"{outdir}/Phyla.png",
    save_params={"w": 1000, "h": 800},
)
tv.show(
    cnidaria,
    legendFun=lambda x: tv.rankLegend(x, "order", ("kingdom",)),
    save_to=f"{outdir}/cnidaria.png",
    save_params={"w": 1000, "h": 800},
)
