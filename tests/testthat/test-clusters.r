library("tidyverse")
# library("go")
available <- BiocManager::available() %>% as.list()
available %>%
  unlist() %>%
  keep(., \(x) str_detect(x, "GO"))
library("glue")

if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  reticulate::use_virtualenv("/home/shannc/.venv/default")
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}

path <- glue("{wd}/tests/test_clusters")
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  sample_name = "A_thaliana",
  combined_results = glue("{path}/data.tsv"),
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo")
)

source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/rrvgo_modified.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))

data <- read_tsv(args$combined_results)
if (!"cluster" %in% colnames(data)) {
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    glue("{path}/embeddings.hdf5"),
    glue("{path}/distances.hdf5")
  )
  py_clusters <- new.env()
  reticulate::source_python(glue("{args$python_source}/clustering.py"), envir = py_clusters)
  dist <- e$cosine
  clusters <- local({
    f <- hclustSk(dist, 0.1, "average", labels_only = FALSE)
    linkage_matrix <- py_clusters$linkageMatrix(f$fitted)
    f$labels
  })

  data <- mergeClusters(clusters, e$metadata, "ProteinId")
}

nested <- data %>%
  group_by(cluster) %>%
  nest() %>%
  mutate(size = map_dbl(data, \(x) nrow(x))) %>%
  arrange(desc(size))

nested[1, ]$data[[1]]$GO_IDs
