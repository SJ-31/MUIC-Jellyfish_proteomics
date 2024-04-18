library(tidyverse)
library(glue)
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}
t <- "tsne"
model <- "prottrans"
args <- list(
  figure_path = glue("{wd}/tests/testthat/output/figures/test_viz_{t}"),
  sample_name = "C_indra",
  compare = TRUE,
  technique = t,
  embedding_path = glue("{wd}/tests/nf-test-out/analysis_subworkflow/Comparison_embeddings/sample_with_comp_embd.hdf5"),
  dist_path = glue("{wd}/tests/nf-test-out/analysis_subworkflow/Comparison_embeddings/sample_with_comp_dist.hdf5"),
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  results_prefix = "testing",
  combined_results = glue("{wd}/results/C_indra_A/1-First_pass/C_indra_all_wcoverage.tsv"),
  uniprot_data = glue("{wd}/data/protein_databases/comparison_taxa/reviewed_all.tsv")
)
source(glue("{wd}/bin/R/analysis/visualize_dr.r"))
main(args)
