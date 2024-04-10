library(fpc)
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
model <- "esm"
MEASURE <- "euclidean"
OUT <- glue("{wd}/tests/testthat/output/clusters_{MEASURE}_{model}")
if (!file.exists(OUT)) {
  dir.create(OUT)
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
  embd_type = "protein",
  sample_embd = glue("{wd}/tests/nf-test-out/C_indra_{model}_embeddings/embeddings.hdf5"),
  sample_embd_dist = glue("{wd}/tests/nf-test-out/C_indra_{model}_embeddings/distances.hdf5")
)

source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
reticulate::use_condaenv(env)

SAMPLE_ONLY <- TRUE
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))

#' BM: Benchmarking functions for clustering methods intended to determine
#' the best parameters, and to compare clustering methods against each
#' other
#' @param dist_t a distance matrix of the embeddings
#'
hierarchichalBM <- function(dist_t, method, height_vec) {
  trees <- list()
  clusters <- hclust(dist_t, method = method)
  for (height in height_vec) {
    cut <- tryCatch(
      expr = cutree(clusters, h = height),
      error = \(cnd) NULL)
    if (!is.null(cut)) {
      trees[[as.character(height)]] <- cut
    }
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

#' Save cluster memberships from a cluster list into a tb, where
#' each column denotes the membership of a given protein under that
#' specific clustering method
saveClusters <- function(cluster_list, method, var, previous_saved) {
  for (cluster in names(cluster_list)) {
    colname <- glue("{method}_{var}_{cluster}")
    dat <- cluster_list[[cluster]]
    tb <- tibble(ProteinId = names(dat),
                 !!colname := dat)
    if (is.null(previous_saved)) {
      previous_saved <- tb
    } else {
      previous_saved <- inner_join(previous_saved, tb, by = join_by(ProteinId))
    }
  }
  return(previous_saved)
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
sample_dist <- sample_protein[[measure]]
ids <- rownames(sample_dist)
sample_dist <- dist(sample_dist)

# Benchmark hierarchichal clustering & dbscan
test_hc <- hierarchichalBM(sample_dist, "average",
                           seq(1, 6, by = 1))
cluster_stats <- getClusterMetrics(sample_dist, test_hc, "hierarchichal_clustering", NULL, "height")
cluster_members <- saveClusters(test_hc, "hc",
                                "height", NULL)

test_hdbscan <- hdbscanBM(sample_dist, seq(5, 30, by = 5))
cluster_stats <- getClusterMetrics(sample_dist, test_hdbscan, "hdbscan", cluster_stats, "min_points")
cluster_members <- saveClusters(test_hc, "hdbscan",
                                "min_points", cluster_members)


# Mask relationships with euclidean distances
# greater than the 3rd quartile for leiden
dist_summary <- summary(sample_dist)
sample_dist[sample_dist > dist_summary[["3rd Qu."]]] <- 0


# Benchmark leiden
reticulate::source_python(glue("{args$python_source}/leiden_clust.py"))
la <- reticulate::import("leidenalg")
py$graph <- createGraph(as.matrix(sample_dist), ids)
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
cluster_members <- saveClusters(test_la, "leiden",
                                "quality_fun", cluster_members)

write_tsv(cluster_stats, glue("{OUT}/cluster_stats.tsv"))
write_tsv(cluster_members, glue("{OUT}/cluster_members.tsv"))
write_tsv(la_metrics, glue("{OUT}/leiden_metrics.tsv"))

