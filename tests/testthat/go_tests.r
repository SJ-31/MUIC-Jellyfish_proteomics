library(tidyverse)
library(glue)
# File for testing out GO analysis
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
  combined_results = glue("{wd}/results/C_indra_A/1-First_pass/C_indra_all_wcoverage.tsv"),
  ontologizer_path = glue("{wd}/tests/nf-test-out/ontologizer/"),
  embedding_path = glue("{wd}/data/reference/go_embedded.npz"),
  dist_path = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"
  ))

## Load samples
source("./bin/R/GO_helpers.r")
source("./bin/R/analysis/metric_functions.r")
source("./bin/R/rrvgo_modified.r")
orgdb_pth <- "./tests/testthat/output/org.Cindrasaksajiae.eg.db"
db_name <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)
rrvgo_path <- "./tests/testthat/output/rrvgo"

# orgdb <- prepOrgDb(orgdb_pth)
# semdata <- lapply(ONTOLOGIES, \(x) {
#   GOSemSim::godata(OrgDb = db_name, ont = x, keytype = "GID")
# }) %>% `names<-`(ONTOLOGIES)

uniprot <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa/reviewed_all.tsv"

sample_name <- "C_indra"


d <- goData(args$combined_results,
            onto_path = args$ontologizer_path)
tb <- read_tsv(args$combined_results)


flattenJoined <- function(vec, split) {
  vec %>%
    discard(is.na) %>%
    lapply(., \(x) str_split_1(x, split)) %>%
    unlist()
}

