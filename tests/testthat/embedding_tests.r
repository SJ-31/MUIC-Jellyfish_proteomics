source("./bin/R/GO_helpers.r")
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}

model <- "a2v"
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  embd_type = "protein",
  sample_name = "C_indra",
  uniprot_data_dir = glue("{wd}//data/protein_databases/comparison_taxa"),
  combined_results = glue("{wd}//results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv"),
  ontologizer_path = glue("{wd}//tests/nf-test-out/ontologizer/"),
  embedding_path = glue("{wd}//tests/nf-test-out/C_indra_{model}_embeddings/embeddings.hdf5"),
  dist_path = glue("{wd}//tests/nf-test-out/C_indra_{model}_embeddings/distances.hdf5"
  ))
source("./bin/R/DR_helpers.r")
source("./bin/R/analysis/prepare_embeddings.r")

embeddingData(args$combined_results,
              args$sample_name,
              args$embedding_path,
              args$dist_path, TRUE)


