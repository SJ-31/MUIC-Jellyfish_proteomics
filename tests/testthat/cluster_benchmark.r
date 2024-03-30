library(clv)
library(fpc)
library(clValid)
library(tidyverse)
library(glue)
library(dbscan)

if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}

args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin/"),
  uniprot_embeddings = glue("{wd}/data/protein_databases/uniprot_embeddings_mean.tsv"),
  uniprot_data_dir = glue("{wd}/data/protein_databases/comparison_taxa"),
  figure_path = glue("{wd}/tests/testthat/output/figures"),
  combined_results = glue("{wd}/results/C_indra/1-First_pass/C_indra_all.tsv"),
  embeddings_path = glue("{wd}/data/reference/go_embedded.npz"),
  ontologizer_path = glue("{wd}/tests/nf-test-out/ontologizer/"),
  sample_name = "C_indra",
  protein_embd_mode = "mean" # One of mean or sum
)

source("./bin/R/GO_helpers.r")
source("./bin/R/DR_helpers.r")
reticulate::use_condaenv(env)

SAMPLE_ONLY <- TRUE
source("./bin/R/analysis/prepare_embeddings.r")

#' BM: Benchmarking functions for clustering methods intended to determine
#' the best parameters, and to compare clustering methods against each
#' other
#' @param dist_t a distance matrix of the embeddings
#'
hierarchichalBM <- function(dist_t, method, height_vec) {
  trees <- list()
  clusters <- hclust(dist_t, method = method)
  for (height in height_vec) {
    trees[[as.character(height)]] <- cutree(clusters, h = height)
  }
  return(trees)
}

hdbscanBM <- function(dist_t, min_points) {
  scans <- list()
  for (m in min_points) {
    hdb <- dbscan::hdbscan(dist_t, minPts = m)
    membership <- hdb$cluster + 1 # fpc stats uses different indices
    names(membership) <- names(hdb$coredist)
    scans[[as.character(m)]] <- membership
  }
  return(scans)
}


#' Compute clustering metrics, appending them to a list
#'
getClusterMetrics <- function(dist_t, cluster_list, method, metrics, var) {
  for (n in names(cluster_list)) {
    stats <- tryCatch(
      expr = fpc::cluster.stats(dist_t, clustering = cluster_list[[n]]),
      error = \(cnd) NULL
    )
    if (!is.null(stats)) {
      current <- discard(stats, \(x) length(x) > 1 || is.null(x))
      current$parameter <- paste0(var, "=", n)
      current$method <- method
      if (is.null(metrics)) {
        metrics <- as.data.frame(current)
      } else {
        metrics <- dplyr::bind_rows(metrics, as.data.frame(current))
      }
    }
  }
  return(metrics)
}


#' Properly format cluster membership from a leidenalg.VertexPartition
#' object i.e. converting it to a named vector
#'
getPartition <- function(py_ptition) {
  vec <- py_ptition$membership + 1
  names(vec) <- py$ptition$graph$vs["name"]
  return(vec)
}


# decrease size for testing
# test <- sample_prot_embd %>% slice(1:500)

sample_dist <- test %>%
  distinct() %>%
  t2Df(., "ProteinId") %>%
  dist()


# Benchmark hierarchichal clustering & dbscan
test_hc <- hierarchichalBM(sample_dist, "average",
                           seq(1, 6, by = 1))
cluster_stats <- getClusterMetrics(sample_dist, test_hc, "hierarchichal_clustering", NULL, "height")
test_hdbscan <- hdbscanBM(sample_dist, seq(5, 30, by = 5))
cluster_stats <- getClusterMetrics(sample_dist, test_hdbscan, "hdbscan", cluster_stats, "min_points")


# Mask relationships with euclidean distances
# greater than the 3rd quartile for leiden
dist_summary <- summary(sample_dist)
sample_dist[sample_dist > dist_summary[["3rd Qu."]]] <- 0


# Benchmark leiden
reticulate::source_python(glue("{args$python_source}/leiden_clust.py"))
la <- reticulate::import("leidenalg")
py$graph <- createGraph(as.matrix(sample_dist), test$ProteinId)
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
}

cluster_stats <- getClusterMetrics(sample_dist, test_la, "leiden", cluster_stats, "quality_fun")

write_tsv(cluster_stats, "./cluster_stats.tsv")
write_tsv(la_metrics, "./leiden_metrics.tsv")

