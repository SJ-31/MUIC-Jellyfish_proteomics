library(fpc)
library(tidyverse)
library(glue)
library(dbscan)

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
    message(glue("hclust at height {height} completed"))
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
    message(glue("hdbscan at min_points {m} completed"))
  }
  return(scans)
}

silhouetteScore <- function(dist_t, cluster_labels) {
  smc <- reticulate::import("sklearn.metrics.cluster")
  result <- smc$silhouette_score(dist_t, cluster_labels, metric = "precomputed")
  return(result)
}

internalMetrics <- function(dist_t, cluster_labels) {
  stats <- list()
  try(stats <- fpc::cluster.stats(dist_t, clustering = cluster_labels, silhouette = FALSE))
  stats$silhouette <- silhouetteScore(dist_t, cluster_labels)
  if (length(stats) == 0) {
    return(NULL)
  }
  return(stats)
}


#' Compute clustering metrics, appending them to a list
#'
getClusterMetrics <- function(dist_t, cluster_labels, method, metrics, var) {
  for (n in names(cluster_labels)) {
    current <- internalMetrics(dist_t, cluster_labels[[n]])
    current$parameter <- paste0(var, "=", n)
    current$method <- method
    if (is.null(metrics)) {
      metrics <- as.data.frame(current)
    } else {
      metrics <- dplyr::bind_rows(metrics, as.data.frame(current))
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
