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
hclust_sk <- function(dist, height, linkage, structured = FALSE, labels_only = TRUE) {
  sc <- reticulate::import("sklearn.cluster")
  if (structured) {
    neighbors <- reticulate::import("sklearn.neighbors")
    knn_graph <- neighbors$kneighbors_graph(dist, as.integer(5), metric = "cosine")
    hc <- sc$AgglomerativeClustering(
      distance_threshold = height,
      n_clusters = NULL,
      linkage = linkage,
      metric = "precomputed",
      connectivity = knn_graph
    )
  } else {
    hc <- sc$AgglomerativeClustering(
      distance_threshold = height,
      n_clusters = NULL,
      linkage = linkage,
      metric = "precomputed"
    )
  }
  if (labels_only) {
    membership <- hc$fit_predict(dist)
    names(membership) <- rownames(dist)
    return(membership + 1)
  }
  fit <- hc$fit(dist)
  labels <- fit$labels_ + 1
  names(labels) <- rownames(dist)
  return(list(fitted = fit, labels = labels))
}


`_hdbscan` <- function(dist, min_points) {
  hdb <- dbscan::hdbscan(dist, minPts = min_points)
  membership <- hdb$cluster + 1 # fpc stats uses different indices
  names(membership) <- names(hdb$coredist)
  return(membership)
}


GRAPH_CREATED <- FALSE
leiden_create_graph <- function(dist) {
  if (!GRAPH_CREATED) {
    GRAPH_CREATED <<- TRUE
    ids <- rownames(dist)
  }
  return(createGraph(as.matrix(dist), ids))
}

LEIDEN_METRICS <- data.frame()
`_leiden` <- function(dist, partition_type) {
  reticulate::source_python(glue("{args$python_source}/clustering.py"))
  if (!GRAPH_CREATED) {
    stop("Need to create graph first!")
  }
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
#' and the value of the parameter being changed
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

sklearn_internal_metrics <- function(dist, cluster_labels, method) {
  smc <- reticulate::import("sklearn.metrics.cluster")
  if (method == "silhouette") {
    return(smc$silhouette_score(dist, cluster_labels, metric = "precomputed"))
  } else if (method == "calinhara") {
    return(smc$calinski_harabasz_score(dist, cluster_labels))
  } else if (method == "davies_bouldin") {
    return(smc$davies_bouldin_score(dist, cluster_labels))
  }
}

show_error <- function(expression, if_failed = NA, label = NULL) {
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


internal_metrics_all <- function(dist_t, cluster_labels) {
  n_unique_labels <- length(unique(cluster_labels))
  if (n_unique_labels == 1) { # don't bother when clustering failed
    return(NULL)
  } else if (n_unique_labels == length(rownames(dist_t))) {
    return(NULL) # Another way the clustering fails is if every
    # element is in a unique cluster
  }
  dunn <- show_error(clValid::dunn(distance = dist_t, clusters = cluster_labels),
    label = "dunn calculation"
  )
  connectivity <- show_error(clValid::connectivity(distance = sample, clusters = dist_matrices),
    label = "connectivity calculation"
  )
  stats <- list(
    silhouette_width = show_error(sklearn_internal_metrics(dist_t, cluster_labels, "silhouette")),
    calinhara = show_error(sklearn_internal_metrics(dist_t, cluster_labels, "calinhara")),
    davies_bouldin = show_error(sklearn_internal_metrics(dist_t, cluster_labels, "davies_bouldin")),
    dunn = dunn,
    connectivity = connectivity
  )
  return(stats)
}

#' Compute clustering metrics, appending them to a list
#' @param cluster_labels: a list of different clusterings, varing by "var"
get_cluster_metrics <- function(dist_t, cluster_labels, method, metrics, var) {
  for (n in names(cluster_labels)) {
    current <- internal_metrics_all(dist_t, cluster_labels[[n]])
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
get_partition <- function(py_ptition) {
  vec <- py_ptition$membership + 1
  names(vec) <- py$ptition$graph$vs["name"]
  return(vec)
}


#' Save cluster memberships from a cluster list into a tb, where
#' each column denotes the membership of a given protein under that
#' specific clustering method
save_clusters <- function(cluster_list, method, var, previous_saved, join_col) {
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


add_cluster_col <- function(args, tb) {
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    args$protrans_embd,
    args$protrans_dist,
  )
  dist <- e$cosine
  clusters <- local({
    labels <- hclust_sk(dist, 0.1, "average")
    tibble(ProteinId = names(labels), cluster = labels)
  })
  tb <- inner_join(tb, clusters)
  return(tb)
}

#' Merge the cluster information in `cluster_labels` with the data in
#' `tb` by `id_col`
merge_clusters <- function(cluster_labels, tb, id_col) {
  assertArg(cluster_labels, \(x) !is.null(names(x)))
  col <- tibble(!!as.symbol(id_col) := names(cluster_labels),
    cluster = as.double(cluster_labels)
  )
  inner_join(tb, col)
}

#' Determine which KEGG pathways strongly overlap with a cluster.
#' @param threshold Overlap threshold
#' @param nested_tb data tb after grouping by then nesting on the column
#' containing the cluster labels
associated_pathways <- function(tb, nested_tb, threshold = 50) {
  kegg_pathways <- groupPathways(tb)
  oneSet <- function(proteins) {
    check <- kegg_pathways$grouped %>%
      lapply(., \(x) {
        if (length(x) < length(proteins)) {
          total_members <- length(x)
        } else {
          total_members <- length(proteins)
        }
        n_intersected <- base::intersect(x, proteins) %>% length()
        percent_in_pathway <- (n_intersected / total_members) * 100
        if (percent_in_pathway > threshold) {
          return(TRUE)
        }
        return(FALSE)
      }) %>%
      unlist()
    kegg_pathways$grouped[check]
  }

  lapply(nested_tb$data, \(x) {
    pathway_names <- oneSet(x$ProteinId) %>% names()
    if (length(pathway_names) == 0) {
      return(NA)
    }
    paste0(pathway_names, collapse = ";")
  }) %>%
    unlist()
}

#' Aggregate identified proteins by specified grouping column
#' @param data A data frame containing the data to be aggregated.
#' @param grouping_col The column by which to group the data.
#' @return A data frame with aggregated metadata, including the size
#' of each group, GO slims, Percolator groups, organisms, header words,
#' GO slim counts, GO IDs, and GO counts.
#' Additionally, KEGG pathways are associated with the groups.
#' @details The function performs the following aggregations:
#'  Groups the data by the specified column.
#'  Computes the size of each group.
#'  Extracts unique GO slims for each group.
#'  Concatenates unique values from the 'Group' column.
#'  Concatenates unique values from the 'organism' column.
#'  Extracts the top 3 words from the 'header' column.
#'  Counts the number of unique GO slims for each group.
#'  Extracts unique GO IDs for each group.
#'  Counts the number of unique GO IDs for each group.
#'  Associates KEGG pathways with the groups.
aggregate_metadata <- function(data, grouping_col) {
  data <- mergeLfq(data, "mean") %>%
    inner_join(., data, by = join_by(ProteinId))
  nested <- data %>%
    group_by(!!as.symbol(grouping_col)) %>%
    nest() %>%
    mutate(
      size = map_dbl(data, \(x) nrow(x)),
      GO_slims = lapply(
        data,
        \(x) get_go_vec(x, go_column = "GO_slims", unique = TRUE)
      ),
      organisms = map_chr(data, \(x) paste0(unique(x$organism), collapse = ";")),
      header_words = map_chr(data, \(x) {
        paste0(header_top_n(x$header, 3), collapse = ";")
      }),
      GO_slim_counts = map_dbl(GO_slims, \(x) length(x)),
      GO_IDs = lapply(
        data,
        \(x) get_go_vec(x, go_column = "GO_IDs", unique = TRUE)
      ),
      GO_counts = map_dbl(GO_IDs, \(x) length(x)),
      GO_category_BP = map_chr(data, \(x) paste0(unique(x$GO_category_BP), collapse = ";")),
      GO_category_CC = map_chr(data, \(x) paste0(unique(x$GO_category_CC), collapse = ";")),
      GO_category_MF = map_chr(data, \(x) paste0(unique(x$GO_category_MF), collapse = ";")),
      mean_coverage = map_dbl(data, \(x) mean(x$pcoverage_nmatch, na.rm = TRUE)),
      median_coverage = map_dbl(data, \(x) median(x$pcoverage_nmatch, na.rm = TRUE)),
      mean_intensity = map_dbl(data, \(x) mean(x$log_intensity, na.rm = TRUE))
    ) %>%
    arrange(desc(size))
  if (grouping_col != "Group") {
    nested <- nested %>% mutate(
      Percolator_groups = map_chr(data, \(x) paste0(unique(x$Group), collapse = ";"))
    )
  }
  nested$KEGG_Pathway <- associated_pathways(data, nested)
  nested
}


REPLACEMENTS <- local({
  original <- c(
    "predicted protein",
    "hypothetical protein",
    "unnamed protein product",
    "light chain",
    "outer membrane protein",
    "heat shock protein",
    "uncharacterized proetin"
  )
  compounded <- map_chr(original, \(x) str_replace_all(x, " ", "_"))
  names(compounded) <- original
  compounded
})


KEYWORDS <- c(
  "hypothetical protein", "ribosomal protein",
  "cytochrome", "histone", "unnamed protein", "toxin",
  "outer membrane protein", "furin", "catenin", "dynein", "tubulin",
  "septin", "uncharacterized protein", "ribosomal", "chaperone", "predicted",
  "kinase", "myosin", "atpase", "zinc finger", "metalloproteinase",
  "microtubule", "nadh-ubiquinone oxidoreductase",
  "ribonucleoprotein", "venom", "ras", "rab", "actin", "proteasome", "elongation factor",
  "scramblase", "protease", "peptidase", "kinesin", "heat shock",
  "adp-ribosylation", "porin", "collagen", "dynamin", "filamin",
  "gelosin", "glutamate", "glutaryl-coa", "glyceraldehyde", "isocitrate"
)

UNWANTED <- c("partial", "fragment", "NA")

name_from_header <- function(header) {
  if (str_detect(header, "\\|")) {
    ex <- header |> str_extract("(tr|sp)\\|[A-Z0-9\\|_]+\\s(.*)OS=.*", group = 2)
  } else {
    ex <- header |> str_extract(".*\\.1(.*)\\[.*", group = 1)
  }
  str_trim(ex) |>
    str_to_lower()
}

#' Any element of 'text_vec' containing a keyword "k" gets mapped to that
#' keyword.
#' @return A quanteda dictionary that can be used with `tokens_lookup`
#' to remove the redundant elements of text_vec
create_dictionary <- function(keywords, text_vec) {
  keyword_list <- setNames(rep(list(c()), length(keywords)), keywords)
  for (word in keywords) {
    no_match <- text_vec |> discard(~ str_detect(.x, word))
    match <- text_vec |> keep(~ str_detect(.x, word))
    keyword_list[[word]] <- match
    text_vec <- no_match
  }
  return(quanteda::dictionary(keyword_list))
}

#' Return the top n words found in a group of FASTA file headers
header_top_n <- function(header_vec, n) {
  headers_only <- map_chr(header_vec, name_from_header) |> discard(is.na)
  dict <- create_dictionary(KEYWORDS, headers_only)
  tokens <- headers_only |>
    quanteda::tokens(
      remove_punct = TRUE, remove_symbols = TRUE,
      remove_numbers = TRUE
    ) |>
    quanteda::tokens_lookup(dict, exclusive = FALSE) |>
    quanteda::tokens_remove(quanteda::stopwords("english")) |>
    quanteda::tokens_remove(UNWANTED)
  freq_matrix <- tokens |>
    quanteda::dfm() |>
    quanteda.textstats::textstat_frequency()
  freq_matrix$feature[1:n]
}


#' Enrich GO IDs in proteins groups using ontologizer
#' @param nested the nested version of `tb`
#' @param tb tbl df containing the columns "ProteinId" and "GO_IDs", used
#' to create the ID -> GO mapping file for ontologizer, as well as the universe
enrich_groups <- function(tb, nested, ontologizer_path, go_path, group_name = "cluster", min_group_size = 10) {
  filter_enriched <- function(ontologizer_df, n_groups, threshold) {
    if (nrow(ontologizer_df) == 0) {
      return(NA)
    }
    tb <- as_tibble(ontologizer_df) |>
      # mutate(
      #   `p.adjusted` = p.adjust(p.adjusted, n = n_groups),
      #   `is.trivial` = case_when(
      #     `p.adjusted` < threshold ~ FALSE,
      #     .default = TRUE
      #   )
      # ) |>
      filter(!str_detect(name, "obsolete") & !is.trivial) |>
      pluck("ID") |>
      paste0(collapse = ";")
  }

  ont <- new.env()
  reticulate::source_python(glue("{args$python_source}/ontologizer_wrapper.py"), envir = ont)
  O <- reticulate_show_error(ont$Ontologizer(tb, ontologizer_path, args$go_path))
  too_small <- nested |> filter(size < min_group_size) # No point when the
  # clusters don't meet the specified threshold
  nested <- nested |> filter(size >= min_group_size)
  groups <- mapply(\(x, y) {
    result <- list()
    result[[glue("{group_name}_{x}")]] <- y$ProteinId
    result
  }, nested[[group_name]], nested$data) |> discard(~ length(.x) == 0)
  params <- list(`-m` = "Bonferroni-Holm")
  enriched <- O$runAll(groups, params)

  enriched <- enriched %>% lapply(., filter_enriched, n_groups = nrow(nested), threshold = 0.05)
  nested$enriched_GO_IDs <- enriched
  too_small$enriched_GO_IDs <- NA
  bind_rows(nested, too_small)
}

#' Save aggregated data to a file
save_aggregated <- function(nested, filename) {
  concatCol <- function(col) {
    lapply(col, \(x) paste0(x, collapse = ";")) |> unlist()
  }
  not_scalars <- colnames(nested)[!sapply(nested, is_scalar_col)]
  nested |>
    select(-data) |>
    mutate(across(contains(not_scalars), concatCol)) |>
    ungroup() |>
    write_tsv(filename)
}
