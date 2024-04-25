library(fpc)
library(tidyverse)
library(glue)
library(dbscan)
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
prepOrgDb(orgdb_pth)

DIST_TYPE <- "a2v"
EMBEDDING_TYPE <- "go"
# LOG
#

d <- goData(args$combined_results)
gos <- d$go_vec$sample

if (DIST_TYPE == "semantic" && EMBEDDING_TYPE == "go") {
  dist <- rrvgo::calculateSimMatrix(gos, orgdb = db_name,
                                    ont = "BP", method = "Wang",
                                    keytype = "GID") %>%
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

# decrease size for testing
#
# Hierarchichal clustering
heights <- seq(min(dist), max(dist), length.out = 12)[-1]
heights <- round(heights[-length(heights)], 2)

test_hc <- benchmarker(dist, `_hclustSkLearn`, as.double(heights), "height", "hclust")
cluster_members <- saveClusters(test_hc, "hc", "height",
                                NULL, JOIN_COL)
write_tsv(cluster_members, glue("{OUTDIR}/cluster_members.tsv"))
cluster_stats <- getClusterMetrics(dist, test_hc, "hierarchichal_clustering", NULL, "height")
write_tsv(cluster_stats, glue("{OUTDIR}/cluster_stats.tsv"))

# hdbscan
test_hdbscan <- benchmarker(dist, `_hdbscan`, seq(5, 30, by = 5), "min_points", "hbscan")
cluster_members <- saveClusters(test_hc, "hdbscan",
                                "min_points", cluster_members,
                                JOIN_COL)
write_tsv(cluster_members, glue("{OUTDIR}/cluster_members.tsv"))
cluster_stats <- getClusterMetrics(dist, test_hdbscan, "hdbscan", cluster_stats, "min_points")
write_tsv(cluster_stats, glue("{OUTDIR}/cluster_stats.tsv"))

# # Mask relationships with euclidean distances
## greater than the 3rd quartile for leiden
dist_summary <- summary(dist)
sample_dist[sample_dist > dist_summary[["3rd Qu."]]] <- 0
#
#
# # Benchmark leiden
reticulate::source_python(glue("{args$python_source}/leiden_clust.py"))

test_la <- benchmarker(dist, `_leiden`,
                       c("Modularity", "RBER", "RB", "CPM", "Surprise"),
                       "partition_type", "leiden")
cluster_stats <- getClusterMetrics(dist, test_la, "leiden", cluster_stats, "quality_fun")
cluster_members <- saveClusters(test_la, "leiden",
                                "quality_fun", cluster_members,
                                JOIN_COL)

write_tsv(cluster_stats, glue("{OUTDIR}/cluster_stats.tsv"))
write_tsv(cluster_members, glue("{OUTDIR}/cluster_members.tsv"))
write_tsv(LEIDEN_METRICS, glue("{OUTDIR}/leiden_metrics.tsv"))
