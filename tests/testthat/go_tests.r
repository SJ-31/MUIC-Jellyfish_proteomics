library(tidyverse)

library(glue)
# File for testing out GO analysis
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  reticulate::use_virtualenv("/home/shannc/.venv/default")
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
  dist_path = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"),
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo")
)

## Load samples
source("./bin/R/GO_helpers.r")
source("./bin/R/analysis/metric_functions.r")
source("./bin/R/rrvgo_modified.r")
orgdb_pth <- "./tests/testthat/output/org.Cindrasaksajiae.eg.db"
DB_NAME <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)
rrvgo_path <- "./tests/testthat/output/rrvgo"

# orgdb <- prepOrgDb(orgdb_pth)
# semdata <- lapply(ONTOLOGIES, \(x) {
#   GOSemSim::godata(OrgDb = db_name, ont = x, keytype = "GID")
# }) %>% `names<-`(ONTOLOGIES)

uniprot <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa/reviewed_all.tsv"

sample_name <- "C_indra"


d <- goData(args$combined_results,
  onto_path = args$ontologizer_path
)
tb <- read_tsv(args$combined_results)


flattenJoined <- function(vec, split) {
  vec %>%
    discard(is.na) %>%
    lapply(., \(x) str_split_1(x, split)) %>%
    unlist()
}

id_with_open <- d$ontologizer$id_with_open_GO
id_with_open <- id_with_open[1:1000]

id_with_open %>%
  map_dbl(., \(x) -log(x)) %>%
  sort()
reduced <- reduceGOList(-log(id_with_open))
reduced$scatter_plots
reduced$reduced_matrix$MF$cluster %>% table()
# reduced$reduced_matrix %>%
#   lapply(., dim) %>%
#   unlist() %>%
