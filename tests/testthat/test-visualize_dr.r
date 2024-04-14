library(tidyverse)
library(glue)
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}
t <- "umap"
model <- "prottrans"
args <- list(
  figure_path = glue("./output/figures/test_viz_{t}"),
  sample_name = "C_indra",
  compare = FALSE,
  technique = t,
  embedding_path = glue("../nf-test-out/C_indra_{model}_embeddings/embeddings.hdf5"),
  dist_path = glue("../nf-test-out/C_indra_{model}_embeddings/distances.hdf5"),
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  results_prefix = "testing",
  combined_results = glue("{wd}/results/C_indra_B/1-First_pass/C_indra_all_wcoverage.tsv"),
  uniprot_data = glue("{wd}/data/protein_databases/comparison_taxa/reviewed_all.tsv")
)
source(glue("{wd}/bin/R/analysis/visualize_dr.r"))
main(args)
