library(tidyverse)
library(glue)
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
  embedding_path = glue("{wd}/data/reference/go_embedded.npz"),
  sample_name = "C_indra",
  uniprot_tsv = glue("{wd}/data/protein_databases/comparison_taxa/reviewed_all.tsv")
)
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))

output <- glue("{wd}/data/protein_databases/uniprot_embeddings/a2v_embeddings.hdf5")

uniprot <- getUniprotData(args$uniprot_tsv)
combineGOEmbd(output, uniprot$go_vec, uniprot$map, NULL)
