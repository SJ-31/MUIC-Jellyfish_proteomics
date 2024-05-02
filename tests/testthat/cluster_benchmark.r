library(fpc)
library(tidyverse)
library(glue)
library(dbscan)
TEST <- FALSE
HIERARCHICHAL_CLUSTERING <- TRUE
HDBSCAN <- TRUE
LEIDEN <- TRUE
# Compare clustering GO terms using semantic distance or cosine distances
# with GO embeddings

if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  setwd("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/testthat")
} else {
  setwd("/mnt/data/shannc/nf/tests/testthat")
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}

args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin/"),
  figure_path = glue("{wd}/tests/testthat/output/figures"),
  combined_results = glue("{wd}/results/C_indra_A/1-First_pass/C_indra_all.tsv"),
  embedding_path = glue("{wd}/tests/nf-test-out/C_indra_a2v_go_embeddings/distances.hdf5"),
  sample_name = "C_indra"
)

source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
orgdb_pth <- "./output/org.Cindrasaksajiae.eg.db"
db_name <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)

DIST_TYPE <- "semantic"
EMBEDDING_TYPE <- "protein"
# LOG
# 2024-05-01 Finished with dist_type = "a2v" and embedding_type = "go"
# 2024-05-01 Finished with dist_type = "semantic" and embedding_type = "go"

d <- goData(args$combined_results)
gos <- d$go_vec$sample

if (DIST_TYPE == "semantic" && EMBEDDING_TYPE == "go") {
  prepOrgDb(orgdb_pth)
  dist <- rrvgo::calculateSimMatrix(gos,
    orgdb = db_name,
    ont = "BP", method = "Wang",
    keytype = "GID"
  ) %>%
    simToDist()
} else if (EMBEDDING_TYPE == "go") {
  dpath <- "../nf-test-out/C_indra_a2v_go_embeddings/distances.hdf5"
  epath <- "../nf-test-out/C_indra_a2v_go_embeddings/embeddings.hdf5"
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    epath,
    dpath
  )
  dist <- e$cosine
} else if (EMBEDDING_TYPE == "protein") {
  dpath <- "../nf-test-out/C_indra_prottrans_embeddings/distances.hdf5"
  epath <- "../nf-test-out/C_indra_prottrans_embeddings/embeddings.hdf5"
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    epath,
    dpath
  )
  dist <- e$cosine
}

if (EMBEDDING_TYPE == "go") {
  OUTDIR <- glue("./output/cluster_go_{DIST_TYPE}")
  JOIN_COL <- "id"
} else {
  OUTDIR <- glue("./output/cluster_prottrans")
  JOIN_COL <- "ProteinId"
}
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR)
}

if (TEST) {
  sample <- rownames(dist) %>% base::sample(size = 500)
  dist <- filterDistMatrix(dist, sample)
}

# decrease size for testing
CLUSTER_STATS <- NULL
CLUSTER_MEMBERS <- NULL
#
# Hierarchichal clustering

if (HIERARCHICHAL_CLUSTERING) {
  heights <- seq(min(dist), max(dist), length.out = 12)[-1]
  heights <- round(heights[-length(heights)], 2)
  test_hc <- benchmarker(dist, `_hclustSkLearn`, as.double(heights), "height", "hclust")
  CLUSTER_MEMBERS <- saveClusters(
    test_hc, "hc", "height",
    CLUSTER_MEMBERS, JOIN_COL
  )
  write_tsv(CLUSTER_MEMBERS, glue("{OUTDIR}/cluster_members.tsv"))
  CLUSTER_STATS <- getClusterMetrics(dist, test_hc, "hierarchichal_clustering", CLUSTER_STATS, "height")
  write_tsv(CLUSTER_STATS, glue("{OUTDIR}/cluster_stats.tsv"))
}

# hdbscan
if (HDBSCAN) {
  test_hdbscan <- benchmarker(dist, `_hdbscan`, seq(5, 30, by = 5), "min_points", "hbscan")
  CLUSTER_MEMBERS <- saveClusters(
    test_hc, "hdbscan",
    "min_points", CLUSTER_MEMBERS,
    JOIN_COL
  )
  write_tsv(CLUSTER_MEMBERS, glue("{OUTDIR}/cluster_members.tsv"))
  CLUSTER_STATS <- getClusterMetrics(dist, test_hdbscan, "hdbscan", CLUSTER_STATS, "min_points")
  write_tsv(CLUSTER_STATS, glue("{OUTDIR}/cluster_stats.tsv"))
}

if (LEIDEN) {
  # # Mask relationships with euclidean distances
  ## greater than the 3rd quartile for leiden (just set them to the maximum distance away)
  dist[dist > quantile(dist, 0.75)] <- max(dist)
  # # Benchmark leiden
  reticulate::source_python(glue("{args$python_source}/leiden_clust.py"))
  py$graph <- leidenCreateGraph(dist)
  test_la <- benchmarker(
    dist, `_leiden`,
    c("Modularity", "RBER", "RB", "CPM", "Surprise"),
    "partition_type", "leiden"
  )
  CLUSTER_STATS <- getClusterMetrics(dist, test_la, "leiden", CLUSTER_STATS, "quality_fun")
  CLUSTER_MEMBERS <- saveClusters(
    test_la, "leiden",
    "quality_fun", CLUSTER_MEMBERS,
    JOIN_COL
  )
  if (!is.null(CLUSTER_STATS)) {
    write_tsv(CLUSTER_STATS, glue("{OUTDIR}/cluster_stats.tsv"))
  }
  write_tsv(CLUSTER_MEMBERS, glue("{OUTDIR}/cluster_members.tsv"))
  write_tsv(LEIDEN_METRICS, glue("{OUTDIR}/leiden_metrics.tsv"))
}
