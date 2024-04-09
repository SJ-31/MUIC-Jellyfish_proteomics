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

# Combine GO terms assigned to each protein in protein_map using
# combine_func
goEmbedding2Prot <- function(protein_map, embedding_tb, combine_func) {
  purrr::map(names(protein_map), ~{
    dplyr::filter(embedding_tb, GO_IDs %in% protein_map[[.x]]) %>%
      reframe(across(where(is.numeric), combine_func)) %>%
      mutate(ProteinId = .x, .before = V1)
  }) %>%
    bind_rows() %>%
    dplyr::filter(!if_any(contains("V"), is.na)) # bug: This shouldn't be necessary
}

# Wrapper for tsne on "data", removing the non-numeric join_col
.tsne <- function(data, join_col) {
  data <- purrr::pluck(data, "embd")
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
.umap <- function(data, join_col, params) {
  data <- purrr::pluck(data, "embd")
  if (!missing(params) && pluck_exists(params, "umap")) {
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


# Run pca on "data" and join scaled results with join_tb
.pca <- function(data, join_col) {
  data <- purrr::pluck(data, "embd")
  results <- list()
  results$prcomp <- data %>%
    dplyr::select(-!!join_col) %>%
    prcomp()
  rownames(results$prcomp$x) <- data[[join_col]]
  results$ve <- pcaVarExplained(results$prcomp)
  return(results)
}


#' Cosine dissimilarity plotted with PCOA
#' @param dist: a matrix with rownames
#' If "dist" is euclidean, performs PCA
.cosinePCoA <- function(data, join_col) {
  dist <- pluck(data, "cosine")
  pcoa <- vegan::wcmdscale(dist, eig = TRUE, add = TRUE)
  rownames(pcoa$points) <- rownames(dist)
  pcoa$pe <- pcoaPerExplained(pcoa)
  return(pcoa)
}

#' Plot dimensionality reduction results
#'
#' @description
#' Generic function for plotting dimensionality
#' reduction results, highlighting on specific column
#' Generates both 3d and 2d plots
plotDr <- function(to_plot, color_col, path, technique, labels) {
  if (technique == "PCA") {
    prefix <- "PC"
  } else {
    prefix <- "V"
  }
  caption_str <- ifelse(labels$caption == "", "",
                        glue("-{gsub(' ', '_', labels$caption)}")
  )
  biplot <- biplotCustom(to_plot,
                         x = glue("{prefix}1"), y = glue("{prefix}2"),
                         colour_column = color_col,
                         labels = labels
  )
  mySaveFig(biplot, glue("{path}/{technique}_biplot-{color_col}{caption_str}.png"))
  v1 <- to_plot[[glue("{prefix}1")]]
  v2 <- to_plot[[glue("{prefix}2")]]
  v3 <- to_plot[[glue("{prefix}3")]]
  plotly::plot_ly(to_plot,
                  x = v1,
                  y = v2,
                  z = v3, color = ~base::get(color_col),
                  type = "scatter3d",
                  marker = list(size = 5)
  ) %>% mySaveFig(., glue("{path}/{technique}_3d-{color_col}{caption_str}.html"))
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
#' Performs the DR technique specified by "technique" on the
#' proteins in the "sample" dataset and "all" dataset combined with
#' other taxa.
#' optionally write output to long-form tsv file ready to plot with ggplot
#'
#' @param dr_data - a list containing the data to apply DR to,
#' (labelled "data"), the tibble to join on for coloring ("tb"),
#' and the columns to color on ("color")
drWrapper <- function(dr_data, join_on, outdir, name, technique, params) {
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  switch(technique,
         "PCA" = { DR <- .pca },
         "PCOA" = { DR <- .cosinePCoA },
         "UMAP" = { DR <- .umap },
         "TSNE" = { DR <- .tsne }
  )
  metadata <- dr_data$data %>% select(join_on, dr_data$color)
  if (!missing(params)) {
    dr_result <- DR(dr_data, join_on, params)
  } else {
    dr_result <- DR(dr_data, join_on)
  }
  if (pluck_exists(dr_result, "prcomp")) {
    data <- dr_result$prcomp[["x"]]
  } else if (any(class(dr_result) == "wcmdscale")) {
    data <- dr_result$points
    colnames(data) <- colnames(data) %>%
      map_chr(., \(x) gsub("Dim", "V", x))
  } else {
    data <- dr_result
  }
  if (technique == "PCA") {
    write_tsv(dr_result$ve, glue("{outdir}/{name}-var_explained.tsv"))
  } else if (technique == "PCOA") {
    write_tsv(dr_result$pe, glue("{outdir}/{name}-p_explained.tsv"))
  }
  to_plot <- data %>%
    as_tibble() %>%
    mutate(!!join_on := metadata[[join_on]]) %>%
    inner_join(., metadata, by = join_by(!!join_on))
  trustworthiness <- trustworthiness(
    column_to_rownames(dr_data$embd, var = join_on), dr_result)

  write_tsv(to_plot, glue("{outdir}/{name}-DR_results.tsv"))
  cat(trustworthiness, file = glue("{outdir}/{name}-trustworthiness.txt"))
  return(list(to_plot = to_plot,
              result = dr_result,
              trustworthiness = trustworthiness))
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
