from context import view_alignments

if "Bio_SDD" in __file__:
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
else:
    wd = "/home/shannc/workflow"

args = {
    "results_file": f"{wd}/results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv",
    "coverage_threshold": 0.9,
    "alignment_file": f"{wd}/results/C_indra/1-First_pass/aligned_peptides.tsv",
    "peptide_map_file": f"{wd}/results/C_indra/1-First_pass/percolator_peptide_map.tsv",
    "outdir": f"{wd}/tests/pytest/output/alignment_viz",
    "mode": "engine_alignment",
}


view_alignments.main(args)
