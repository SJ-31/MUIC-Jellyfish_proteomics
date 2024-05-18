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
  protrans_dist = glue("{wd}/tests/nf-test-out/C_indra_prottrans_embeddings/distances.hdf5"),
  protrans_embd = glue("{wd}/tests/nf-test-out/C_indra_prottrans_embeddings/embeddings.hdf5"),
  dist_path = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"),
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo")
)

## Load samples
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/rrvgo_modified.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
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
e <- embeddingData(
  args$combined_results,
  args$sample_name,
  args$protrans_embd,
  args$protrans_dist,
)

id_with_open <- d$ontologizer$id_with_open_GO

dist <- e$cosine


OUTDIR <- "/home/shannc/Downloads/thesis_testzone/"
py_clusters <- new.env()
reticulate::source_python(glue("{args$python_source}/clustering.py"), envir = py_clusters)
clusters <- local({
  f <- hclustSk(dist, 0.1, "average", labels_only = FALSE)
  linkage_matrix <- py_clusters$linkageMatrix(f$fitted)
  py_clusters$saveDendogram(linkage_matrix,
    glue("{OUTDIR}/dendogram.png"),
    cutoff = 0.1
  )
  f$labels
})

e$metadata <- mergeClusters(clusters, e$metadata, "ProteinId")

e$metadata %>% colnames()
nested <- e$metadata %>%
  group_by(cluster) %>%
  nest() %>%
  mutate(
    size = map_dbl(data, \(x) nrow(x)),
    GO_slims = lapply(
      data,
      \(x) goVector(x, go_column = "GO_slims", unique = TRUE)
    ),
    GO_slim_counts = map_dbl(GO_slims, \(x) length(x)),
    GO_IDs = lapply(
      data,
      \(x) goVector(x, go_column = "GO_IDs", unique = TRUE)
    ),
    GO_counts = map_dbl(GO_IDs, \(x) length(x))
  ) %>%
  arrange(desc(size))
slims <- lapply(nested$GO_slims, idsIntoOntology) %>% purrr::reduce(., mergeLists)
slim_tb <- nested %>%
  select(cluster, size) %>%
  ungroup() %>%
  mutate(
    cluster = as.vector(cluster),
    slim_CC = slims$CC,
    slim_BP = slims$BP,
    slim_MF = slims$MF
  )

# TODO: Summarize the go terms in each cluster using text-mining techniques
