source("./bin/R/GO_helpers.r")
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  embd_type = "protein",
  sample_name = "C_indra",
  uniprot_data_dir = glue("{wd}/data/protein_databases/comparison_taxa"),
  combined_results = glue("{wd}/results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv"),
  ontologizer_path = glue("{wd}/tests/nf-test-out/ontologizer/"),
  sample_embd = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/embeddings.hdf5"),
  sample_embd_dist = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"
  ))
source("./bin/R/DR_helpers.r")
testthat <- glue("{wd}//tests/testthat/output")

source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
