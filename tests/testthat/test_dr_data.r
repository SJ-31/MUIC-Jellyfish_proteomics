library(tidyverse)
library(glue)

if (EMBEDDING_TYPE == "protein") {
  dpath <- "../nf-test-out/C_indra_a2v_embeddings/distances.hdf5"
  epath <- "../nf-test-out/C_indra_a2v_embeddings/embeddings.hdf5"
  fpath <- "./output/go_dr_protein_level"
} else {
  dpath <- "../nf-test-out/C_indra_a2v_go_embeddings/distances.hdf5"
  epath <- "../nf-test-out/C_indra_a2v_go_embeddings/embeddings.hdf5"
  fpath <- "./output/go_dr_go_level"
}
go_file <- glue("{DATA}/reference/go.obo")
go_slim_file <- glue("{DATA}/reference/goslim_generic.obo")

args <- list(
  combined_results = "../../results/C_indra_A/1-First_pass/C_indra_all.tsv",
  sample_name = "C_indra",
  ontologizer_path = "../nf-test-out/ontologizer/",
  dist_path = dpath,
  embedding_path = epath,
  python_source = "../../bin",
  r_source = "../../bin/R",
  figure_path = fpath
)


techniques <- c("umap", "pcoa", "tsne")
orgdb_pth <- "./output/org.Cindrasaksajiae.eg.db"
db_name <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
prepOrgDb(orgdb_pth)

d <- goData(args$combined_results,
  onto_path = args$ontologizer_path
)

e <- embeddingData(
  args$combined_results,
  args$sample_name,
  args$embedding_path,
  args$dist_path,
)
