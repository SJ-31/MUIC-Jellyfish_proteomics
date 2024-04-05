# Helper functions for analyzing GO embeddings and performing DR

library(tidyverse)
library(glue)
library(umap)
library(tsne)
library(plotly)
library(reticulate)
source(glue::glue("{args$r_source}/rrvgo_modified.r"))

#' Wrapper for sklearn.manifold's trustworthiness function
#'
#' @description
#' X{array-like, sparse matrix} of shape (n_samples, n_features) or
# (n_samples, n_samples) Contains a sample per row.
# X_embedded{array-like, sparse matrix} of shape (n_samples, n_components) Embedding of the training data in low-dimensional space.
# n_neighbors int, default=5
# The number of neighbors that will be considered. Should be fewer than n_samples / 2 to ensure the trustworthiness to lies within [0, 1], as mentioned in [1]. An error will be raised otherwise.
# Considers the first three components by default
trustworthiness <- function(X, X_embedded) {
  if (pluck_exists(X_embedded, "prcomp")) {
    X_embedded <- X_embedded$prcomp$x
  } else if (any(class(X_embedded) == "wcmdscale")) {
    X_embedded <- X_embedded$points
  }
  X_embedded <- X_embedded[, c(1, 2)]
  mf <- import("sklearn.manifold")
  return(mf$trustworthiness(X, X_embedded, n_neighbors = 5L))
}


#' The percentage of variance explained from PCoA
#'
pcoaPerExplained <- function(wcmdscale_obj) {
  pe <- 100 * (wcmdscale_obj$eig / sum(wcmdscale_obj$eig))
  return(tibble(
    PCo = seq(pe),
    Percent_explained = pe,
    Cum_sum_percent_explained = cumsum(pe)
  ))
}


pcaVarExplained <- function(prcomp_obj) {
  vars <- (prcomp_obj$sdev^2 / sum(prcomp_obj$sdev^2))
  return(tibble(
    PC = seq(vars),
    Variance_explained = vars * 100,
    Cum_sum_variance_explained = cumsum(vars * 100)
  ))
}

biplotCustom <- function(ordination_tb, colour_column, x, y, palette, labels) {
  # Convenience function for plotting ordination results
  if (missing(x) || missing(y)) {
    y <- "PC1"
    x <- "PC2"
  }
  if (missing(palette)) {
    p <- "ggprism::colors" # Has 20 colors
  } else {
    p <- palette
  }
  plotted <- ggplot(ordination_tb, aes(
    x = .data[[x]], y = .data[[y]],
    colour = .data[[colour_column]]
  )) +
    geom_point(size = 1, stroke = 1) +
    scale_color_paletteer_d(p) +
    theme_bw() +
    labs(
      title = labels$title,
      caption = labels$caption
    )
  return(plotted)
}

goEmbedding2Prot <- function(protein_map, embedding_tb, combine_func) {
  # Combine GO terms assigned to each protein in protein_map using
  # combine_func
  purrr::map(names(protein_map), ~{
    dplyr::filter(embedding_tb, GO_IDs %in% protein_map[[.x]]) %>%
      reframe(across(where(is.numeric), combine_func)) %>%
      mutate(ProteinId = .x, .before = V1)
  }) %>%
    bind_rows() %>%
    dplyr::filter(!if_any(contains("V"), is.na)) # bug: This shouldn't be necessary
}

# Wrapper for tsne on "data", removing the non-numeric join_col
tsneAndJoin <- function(data, join_col) {
  data_only <- data %>% dplyr::select(-!!join_col)
  tsne <- tsne(data_only, k = 3)
  rownames(tsne) <- data[[join_col]]
  return(tsne)
}

#' Wrapper for Scikit learn's more optimized tsne
tsneSkAndJoin <- function(data, join_col, params) {
  converted <- data %>% column_to_rownames(var = join_col)
  mf <- import("sklearn.manifold")
  if (missing(params)) {
    tsne <- mf$TSNE(n_components = 3L)
  } else {
    tsne <- mf$TSNE(
      n_components = 3L,
      perplexity = params$perplexity,
      metric = params$metric
    )
  }
  result <- tsne$fit_transform(converted)
  rownames(result) <- data[[join_col]]
  return(result)
}

# Perform umap on "data",
umapAndJoin <- function(data, join_col, params) {
  if (!missing(params) && pluck_exists(params, "umap")) {
    print(params)
    umap_params <- params$umap
  } else {
    umap_params <- umap.defaults
    umap_params$n_neighbors <- 10
    umap_params$n_components <- 3
    umap_params$min_dist <- 0.25
  }
  data_only <- data %>% dplyr::select(-!!join_col)
  umap <- umap(data_only, config = umap_params)
  rownames(umap$layout) <- data[[join_col]]
  return(umap$layout)
}

#' Proof of data reconstruction by PCA
#'
#' @description
#' Attempt to reconstruct original data from the pcs
pcaReconstruct <- function(prcomp_obj, num_pcs) {
  if (!missing(num_pcs)) {
    xhat <- prcomp_obj$x[, 1:num_pcs] %*% t(prcomp_obj$rotation[, 1:num_pcs])
  } else {
    xhat <- prcomp_obj$x %*% t(prcomp_obj$rotation)
  }
  if (prcomp_obj$scale) {
    xhat <- scale(xhat, center = FALSE, scale = 1 / prcomp_obj$scale)
  }
  xhat <- scale(xhat, center = -prcomp_obj$center, scale = FALSE)
  return(xhat)
}

# Run pca on "data" and join scaled results with join_tb
pcaAndJoin <- function(data, join_col) {
  results <- list()
  results$prcomp <- data %>%
    dplyr::select(-!!join_col) %>%
    prcomp()
  rownames(results$prcomp$x) <- data[[join_col]]
  results$ve <- pcaVarExplained(results$prcomp)
  return(results)
}

euclideanDistance <- function(v1, v2) {
  return(sqrt(sum((v1 - v2)^2)))
}


#' Cosine dissimilarity plotted with PCOA
#'
cosinePcoaAndJoin <- function(data, join_col) {
  dist <- cosineDissimilarity(data, join_col)
  pcoa <- vegan::wcmdscale(dist, eig = TRUE, add = TRUE)
  rownames(pcoa$points) <- data[[join_col]]
  pcoa$pe <- pcoaPerExplained(pcoa)
  return(pcoa)
}

#' Plot dimensionality reduction results
#'
#' @description
#' Generic function for plotting dimensionality
#' reduction results, highlighting on specific column
#' Generates both 3d and 2d plots
#' Currently works with results from PCA or umap
plotDr <- function(dr_result, join_tb, join_col, color_col, path, name, labels) {
  caption_str <- ifelse(labels$caption == "", "",
                        glue("-{gsub(' ', '_', labels$caption)}")
  )
  prefix <- "V"
  if (pluck_exists(dr_result, "prcomp")) {
    prefix <- "PC"
    dr_result <- dr_result$prcomp[["x"]]
  } else if (any(class(dr_result) == "wcmdscale")) {
    dr_result <- dr_result$points
    colnames(dr_result) <- colnames(dr_result) %>%
      map_chr(., \(x) gsub("Dim", "V", x))
  }
  to_plot <- dr_result %>%
    as_tibble() %>%
    mutate(!!join_col := join_tb[[join_col]]) %>%
    inner_join(., join_tb, by = join_by(!!join_col))
  biplot <- biplotCustom(to_plot,
                         x = glue("{prefix}1"), y = glue("{prefix}2"),
                         colour_column = color_col,
                         labels = labels
  )
  mySaveFig(biplot, glue("{path}/{name}_biplot-{color_col}{caption_str}.png"))
  v1 <- to_plot[[glue("{prefix}1")]]
  v2 <- to_plot[[glue("{prefix}2")]]
  v3 <- to_plot[[glue("{prefix}3")]]
  plotly::plot_ly(to_plot,
                  x = v1,
                  y = v2,
                  z = v3, color = ~base::get(color_col),
                  type = "scatter3d",
                  marker = list(size = 5)
  ) %>% mySaveFig(., glue("{path}/{name}_3d-{color_col}{caption_str}.html"))
  return(biplot)
}

labelGen <- function(analysis_name, sample_name, caption) {
  return(list(
    title = glue("{analysis_name} of {sample_name}"),
    caption = caption
  ))
}

#' Generic dimensionality reduction function
#'
#' @description
#' Performs the DR technique specified by "dR" on the
#' proteins in the "sample" dataset and "all" dataset combined with
#' other taxa
#'
#' @param dr_data - a list containing the data to apply DR to,
#' (labelled "data"), the tibble to join on for coloring ("tb"),
#' and the columns to color on ("color")
#' Note data is a named vector, which will be iterated over
completeDR <- function(dr_data, fig_dir, join_on, prefix, dR, params, title, label) {
  data <- purrr::pluck(dr_data, "data")
  if (!missing(params)) {
    print("parameters passed")
    dr_result <- dR(data, join_on, params)
  } else {
    dr_result <- dR(data, join_on)
  }
  if (missing(label)) {
    label <- labelGen(title, glue("{prefix}"), caption = "")
  }
  dr_data$dr <- dr_result
  if (!is.null(dr_data$embd_type)) {
    path_name <- paste0(fig_dir)
  } else {
    path_name <- paste0(fig_dir, "_", dr_data$embd_type)
  }
  dr_data$biplots <- purrr::map(dr_data$color, \(x) {
    plotDr(dr_result,
           color_col = x,
           labels = label,
           join_tb = dr_data$tb,
           join_col = join_on,
           path = path_name,
           name = paste0(prefix, "_", x)
    )
  })
  dr_data$trustworthiness <- trustworthiness(column_to_rownames(data, var = join_on), dr_result)
  cat(
    dr_data$trustworthiness,
    file = glue("{path_name}/{prefix}_trustworthiness.txt")
  )
  return(dr_data)
}

#' Compute the distance of all rows in a matrix against one vector
#'
#' @description
#' Uses the function specified by "dist_func"
distRows <- function(query, matrix, dist_func) {
  distances <- map(rownames(matrix), \(x) dist_func(matrix[x,], query)) %>%
    `names<-`(rownames(matrix)) %>%
    unlist()
  return(distances)
}

#' Find k nearest proteins in a reduced dimension space, receiving both queries and the points as a named matrix
#'
#' @description
#' Currently supports PCA and PCoA
nearestInDim <- function(query, k, data, dist_func, id_col) {
  if (any(class(data) == "prcomp")) {
    points <- data$x
  } else if (any(class(data) == "wcmdscale")) {
    points <- data$points
  } else {
    points <- data
  }
  if (any(class(data) == "tbl") && !missing(id_col)) {
    points <- column_to_rownames(data, var = id_col)
  }
  result_list <- list()
  for (q in seq_len(nrow(query))) {
    result <- distRows(query[q,], points, dist_func)
    result <- result %>%
      discard_at(\(x) x %in% rownames(query)) %>%
      sort() %>%
      head(n = k)
    result_list[[rownames(query)[q]]] <- result
  }
  result_list <- result_list %>%
    unname() %>%
    unlist() %>%
    uniqueNames()
  return(result_list)
}
