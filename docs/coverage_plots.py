import polars as pl
from matplotlib.figure import Figure
import pandas as pd
import polars.selectors as cs
import sys
from rapidfuzz import distance as di
from pathlib import Path

if "Bio_SDD" in str(Path().absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
else:
    wd = "/home/shannc/workflow"


outdir = f"{wd}/docs/figures/alignments"
python_source = f"{wd}/bin"
sys.path.append(python_source)
from view_alignments import PeptideViz, COLOR_SCHEME
import helpers as hh
import trace_alignments as ta

chosen_pass = "2-Second_pass"
prefixes = {
    "ND": "ND_C_indra",
    "default": "C_indra",
    "Calibrated": "C_indra.calibrated",
}
results = f"{wd}/results"
chosen_path = f"{results}/{prefixes['Calibrated']}/{chosen_pass}"

args = {
    "alignment_file": f"{chosen_path}/aligned_peptides.tsv",
    "peptide_map_file": f"{chosen_path}/percolator_peptide_map_all.tsv",
    "results_file": f"{chosen_path}/{prefixes['Calibrated']}_all_wcoverage.tsv",
}


T = ta.AlignmentTracer(args["alignment_file"], args["peptide_map_file"])
T.get_id2metadata(args["results_file"])

wanted_ids = {
    "cftx2": ("P1045471", 100),  # Toxin CfTX-2
    "ompa": ("P215233", 50),  # OmpA Serratia fonticola
    "ompa_rah": ("P43540", 100),  # OmpA Rahnella contaminans
    "hyla_crisp": ("P5899", 100),  # CRISP from frog
}


for k, v in wanted_ids.items():
    fig: Figure = T.plot_engines_alignment(v[0], wrap=v[1])
    fig.savefig(f"{outdir}/{k}.png", dpi=200, bbox_inches="tight")
