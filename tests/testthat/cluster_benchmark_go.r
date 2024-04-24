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

DIST_TYPE <- "semantic"

d <- goData(args$combined_results)
gos <- d$go_vec$sample

if (DIST_TYPE == "semantic") {
  dist <- rrvgo::calculateSimMatrix(gos, orgdb = db_name,
                                    ont = "MF", method = "Wang",
                                    keytype = "GID")
} else {
  dpath <- "../nf-test-out/C_indra_a2v_go_embeddings/distances.hdf5"
  epath <- "../nf-test-out/C_indra_a2v_go_embeddings/embeddings.hdf5"
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    epath,
    dpath
  )
  dist <- e$cosine
}
OUTDIR <- "./output/cluster_go"

# decrease size for testing
#
# Benchmark hierarchichal clustering & dbscan
test_hc <- hierarchichalBM(dist, "average",
                           seq(1, 6, by = 1))
cluster_members <- saveClusters(test_hc, "hc",
                                "height", NULL)
write_tsv(cluster_members, glue("{OUTDIR}/cluster_members.tsv"))
cluster_stats <- getClusterMetrics(dist, test_hc, "hierarchichal_clustering", NULL, "height")
write_tsv(cluster_stats, glue("{OUTDIR}/cluster_stats.tsv"))

test_hdbscan <- hdbscanBM(dist, seq(5, 30, by = 5))
cluster_members <- saveClusters(test_hc, "hdbscan",
                                "min_points", cluster_members)
write_tsv(cluster_members, glue("{OUTDIR}/cluster_members.tsv"))
cluster_stats <- getClusterMetrics(dist, test_hdbscan, "hdbscan", cluster_stats, "min_points")
write_tsv(cluster_stats, glue("{OUTDIR}/cluster_stats.tsv"))
#
# # Mask relationships with euclidean distances
# # greater than the 3rd quartile for leiden
dist_summary <- summary(dist)
sample_dist[sample_dist > dist_summary[["3rd Qu."]]] <- 0
#
#
# # Benchmark leiden
reticulate::source_python(glue("{args$python_source}/leiden_clust.py"))
la <- reticulate::import("leidenalg")
ids <- rownames(dist)
py$graph <- createGraph(as.matrix(dist), ids)
la_quality_funs <- list("Modularity" = la$ModularityVertexPartition,
                        "RBER" = la$RBERVertexPartition,
                        "RB" = la$RBConfigurationVertexPartition,
                        "CPM" = la$CPMVertexPartition,
                        "Surprise" = la$SurpriseVertexPartition)
la_metrics <- data.frame()
test_la <- list()
for (type in names(la_quality_funs)) {
  py$ptition <- partitionOptimise(py$graph,
                                  la_quality_funs[[type]])
  row <- as.data.frame(partitionMetrics(py$ptition))
  row$type <- type
  la_metrics <- dplyr::bind_rows(la_metrics, row)
  test_la[[type]] <- getPartition(py$ptition)
  message(glue("leiden with partition type {type} completed"))
}
cluster_stats <- getClusterMetrics(dist, test_la, "leiden", cluster_stats, "quality_fun")
cluster_members <- saveClusters(test_la, "leiden",
                                "quality_fun", cluster_members)

write_tsv(cluster_stats, glue("{OUTDIR}/cluster_stats.tsv"))
write_tsv(cluster_members, glue("{OUTDIR}/cluster_members.tsv"))
write_tsv(la_metrics, glue("{OUTDIR}/leiden_metrics.tsv"))
