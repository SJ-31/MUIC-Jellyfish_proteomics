library(fpc)
library(clValid)
library(tidyverse)
library(glue)
library(dbscan)

#' BM: Benchmarking functions for clustering methods intended to determine
#' the best parameters, and to compare clustering methods against each
#' other
#' @param dist_t a distance matrix of the embeddings
#'
`_hclustSkLearn` <- function(dist, height) {
  sc <- reticulate::import("sklearn.cluster")
  hc <- sc$AgglomerativeClustering(
    distance_threshold = height,
    n_clusters = NULL,
    linkage = "average",
    metric = "precomputed"
  )
  membership <- hc$fit_predict(dist)
  names(membership) <- rownames(dist)
  return(membership + 1)
}



`_hclust` <- function(dist, height) {
  clusters <- hclust(dist, method = "average")
  cut <- NULL
  try(cut <- cutree(clusters, h = height))
  return(cut)
}

`_hdbscan` <- function(dist, min_points) {
  hdb <- dbscan::hdbscan(dist, minPts = min_points)
  membership <- hdb$cluster + 1 # fpc stats uses different indices
  names(membership) <- names(hdb$coredist)
  return(membership)
}


GRAPH_CREATED <- FALSE
leidenCreateGraph <- function(dist) {
  if (!GRAPH_CREATED) {
    GRAPH_CREATED <<- TRUE
    ids <- rownames(dist)
  }
  return(createGraph(as.matrix(dist), ids))
}

LEIDEN_METRICS <- data.frame()
`_leiden` <- function(dist, partition_type) {
  if (!GRAPH_CREATED) {
    stop("Need to create graph first!")
  }
  reticulate::source_python(glue("{args$python_source}/leiden_clust.py"))
  la <- reticulate::import("leidenalg")
  la_quality_funs <- list(
    "Modularity" = la$ModularityVertexPartition,
    "RBER" = la$RBERVertexPartition,
    "RB" = la$RBConfigurationVertexPartition,
    "CPM" = la$CPMVertexPartition,
    "Surprise" = la$SurpriseVertexPartition
  )
  py$ptition <- partitionOptimise(py$graph, la_quality_funs[[partition_type]])
  row <- as.data.frame(partitionMetrics(py$ptition))
  row$type <- partition_type
  LEIDEN_METRICS <<- dplyr::bind_rows(LEIDEN_METRICS, row)
  return(getPartition(py$ptition))
}

#' Generic cluster benchmarking function
#'
#' @description
#' @param cluster_fun A function that returns the cluster membership of the
#' elements in dist. This should be a wrapper function on a clustering method
#' that accepts two arguments only: the distance matrix for the clustering
#' and the value of the parameter beign changed
#' @param_vector Vector containing the parameters that will be varied
benchmarker <- function(dist, clusterFun, param_vector, param_name, method_name) {
  clusterings <- list()
  for (p in param_vector) {
    clusters <- clusterFun(dist, p)
    message(glue("{method_name} at {param_name} {p} completed"))
    clusterings[[as.character(p)]] <- clusters
  }
  return(clusterings)
}

sklearnInternalMetrics <- function(dist, cluster_labels, method) {
  smc <- reticulate::import("sklearn.metrics.cluster")
  if (method == "silhouette") {
    return(smc$silhouette_score(dist, cluster_labels, metric = "precomputed"))
  } else if (method == "calinhara") {
    return(smc$calinski_harabasz_score(dist, cluster_labels))
  } else if (method == "davies_bouldin") {
    return(smc$davies_bouldin_score(dist, cluster_labels))
  }
}

showError <- function(expression, if_failed = NA, label = NULL) {
  tryCatch(expr = expression, error = \(cnd) {
    if (!is.null(label)) {
      warn <- glue("Ignored error: {label}")
    } else {
      warn <- "Ignored error"
    }
    message(glue("------\n{warn}\n\tCall: "))
    message(conditionCall(cnd))
    message("------")
    return(if_failed)
  })
}

showError(log("3"))

internalMetricsAll <- function(dist_t, cluster_labels) {
  n_unique_labels <- length(unique(cluster_labels))
  if (n_unique_labels == 1) { # don't bother when clustering failed
    return(NULL)
  } else if (n_unique_labels == length(rownames(dist_t))) {
    return(NULL) # Another way the clustering fails is if every
    # element is in a unique cluster
  }
  dunn <- showError(clValid::dunn(distance = dist_t, clusters = cluster_labels),
    label = "dunn calculation"
  )
  connectivity <- showError(clValid::connectivity(distance = sample, clusters = clusters),
    label = "connectivity calculation"
  )
  stats <- list(
    silhouette_width = sklearnInternalMetrics(dist_t, cluster_labels, "silhouette"),
    calinhara = sklearnInternalMetrics(dist_t, cluster_labels, "calinhara"),
    davies_bouldin = sklearnInternalMetrics(dist_t, cluster_labels, "davies_bouldin"),
    dunn = dunn,
    connectivity = connectivity
  )
  return(stats)
}

#' Compute clustering metrics, appending them to a list
#'
getClusterMetrics <- function(dist_t, cluster_labels, method, metrics, var) {
  for (n in names(cluster_labels)) {
    current <- internalMetricsAll(dist_t, cluster_labels[[n]])
    if (is.null(current)) {
      next
    }
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
saveClusters <- function(cluster_list, method, var, previous_saved, join_col) {
  for (cluster in names(cluster_list)) {
    colname <- glue("{method}_{var}_{cluster}")
    tb <- tibble(
      !!join_col := names(cluster_list[[cluster]]),
      !!colname := as.double(cluster_list[[cluster]])
    )
    if (is.null(previous_saved)) {
      previous_saved <- tb
    } else {
      previous_saved <- inner_join(previous_saved, tb, by = join_by(!!join_col))
    }
  }
  return(previous_saved)
}
