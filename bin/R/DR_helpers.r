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
trustworthiness <- function(X, X_embedded, precomputed = FALSE) {
  n_neighbors <- 5L
  if (pluck_exists(X_embedded, "prcomp")) {
    X_embedded <- X_embedded$prcomp$x
  } else if (any(class(X_embedded) == "wcmdscale")) {
    X_embedded <- X_embedded$points
  }
  X_embedded <- X_embedded[, c(1, 2)]
  mf <- import("sklearn.manifold")
  if (precomputed) {
    return(mf$trustworthiness(X, X_embedded,
      n_neighbors = n_neighbors,
      metric = "precomputed"
    ))
  }
  return(mf$trustworthiness(X, X_embedded, n_neighbors = n_neighbors))
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
  tryCatch(
    expr = {
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
    },
    error = \(cnd) {
      message("Error in plotting")
      print(ordination_tb)
      stop()
    }
  )
  return(plotted)
}


PARAMS <- list(
  tsne = list(n_components = 2, perplexity = 40),
  umap = list(
    n_components = 2,
    n_neighbors = 15,
    min_dist = 1
  )
)
PARAMS$metric <- "cosine"

#' Dimensionality reduction functions for use with generic function
#' With the exception of pca, each function expects `data` to be a list
#' containing a pre-computed distance matrix
#' (either cosine or euclidean)
#' Which measure to use is specified with global variable PARAMS$metric

`_tsne` <- function(data) {
  mf <- reticulate::import("sklearn.manifold")
  if (pluck_exists(data, "cosine")) {
    dist <- purrr::pluck(data, "cosine")
  } else {
    dist <- purrr::pluck(data, "dist")
  }
  tsne <- mf$TSNE(
    n_components = as.integer(PARAMS$tsne$n_components),
    perplexity = as.integer(PARAMS$tsne$perplexity),
    metric = "precomputed",
    init = "random"
  )
  result <- tsne$fit_transform(dist) %>% as.data.frame()
  rownames(result) <- rownames(dist)
  return(result)
}

`_umap` <- function(data) {
  if (pluck_exists(data, "cosine")) {
    dist <- purrr::pluck(data, "cosine")
  } else {
    dist <- purrr::pluck(data, "dist")
  }
  umap <- reticulate::import("umap")
  fit <- umap$UMAP(
    metric = "precomputed",
    n_components = as.integer(PARAMS$umap$n_components),
    n_neighbors = as.integer(PARAMS$umap$n_neighbors),
    min_dist = PARAMS$umap$min_dist
  )
  tryCatch(
    expr = {
      result <- fit$fit_transform(dist) %>% as.data.frame()
    },
    error = \(cnd) {
      message("UMAP fit failed, inspect problematic df as FAILED")
      FAILED <<- dist
      stop()
    }
  )
  rownames(result) <- rownames(dist)
  return(result)
}

simToDist <- function(matrix) {
  max_sim <- max(matrix)
  matrix <- max_sim - matrix
  return(matrix)
}

filterDistMatrix <- function(matrix, labels_to_keep) {
  colnames(matrix) <- rownames(matrix)
  matrix <- matrix[
    rownames(matrix) %in% labels_to_keep,
    colnames(matrix) %in% labels_to_keep
  ]
  return(matrix)
}

`_pca` <- function(data) {
  data <- purrr::pluck(data, "embd")
  results <- list()
  results$prcomp <- prcomp(data)
  rownames(results$prcomp$x) <- rownames(data)
  results$ve <- pcaVarExplained(results$prcomp)
  return(results)
}

`_pcoa` <- function(data) {
  if (pluck_exists(data, "cosine")) {
    dist <- purrr::pluck(data, "cosine")
  } else {
    dist <- purrr::pluck(data, "dist")
  }
  pcoa <- vegan::wcmdscale(dist, eig = TRUE, add = TRUE)
  rownames(pcoa$points) <- rownames(dist)
  pcoa$pe <- pcoaPerExplained(pcoa)
  return(pcoa)
}

#' Generic dimensionality reduction function
#'
#' @description
#' Performs the DR technique specified by "technique" on the
#' proteins in the "sample" dataset and "all" dataset combined with
#' other taxa.
#' Each function called by technique must either return
#' Optionally write output to long-form tsv file ready to plot with ggplot
#'
#' @param dr_data - a list containing the data to apply DR to,
#' (labelled "data"), the tibble to join on for coloring ("tb"),
#' and the columns to color on ("color")
drWrapper <- function(dr_data, join_on, outdir, name, technique) {
  if (!file.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  switch(technique,
    "pca" = {
      DR <- `_pca`
    },
    "pcoa" = {
      DR <- `_pcoa`
    },
    "umap" = {
      DR <- `_umap`
    },
    "tsne" = {
      DR <- `_tsne`
    }
  )
  metadata <- dr_data$metadata %>% dplyr::select(c(join_on, dr_data$color))
  dr_result <- DR(dr_data)
  if (technique == "pca") {
    reduced <- dr_result$prcomp[["x"]]
    write_tsv(dr_result$ve, glue("{outdir}/{name}-var_explained.tsv"))
  } else if (technique == "pcoa") {
    reduced <- dr_result$points
    colnames(reduced) <- map_chr(
      colnames(reduced),
      \(x) gsub("Dim", "V", x)
    )
    write_tsv(dr_result$pe, glue("{outdir}/{name}-p_explained.tsv"))
  } else {
    reduced <- dr_result
  }
  to_plot <- reduced %>%
    as_tibble() %>%
    mutate(!!join_on := metadata[[join_on]]) %>%
    inner_join(., metadata, by = join_by(!!join_on))

  if (purrr::pluck_exists(dr_data, "embd")) {
    trustworthiness <- trustworthiness(dr_data$embd, dr_result)
  } else {
    trustworthiness <- trustworthiness(dr_data$dist, dr_result,
      precomputed = TRUE
    )
  }
  cat(trustworthiness, file = glue("{outdir}/{name}-trustworthiness.txt"))

  write_tsv(to_plot, glue("{outdir}/{name}-DR_results.tsv"))
  return(list(
    to_plot = to_plot,
    result = dr_result,
    trustworthiness = trustworthiness
  ))
}



#' Compute the distance of all rows in a matrix against one vector
#'
#' @description
#' Uses the function specified by "dist_func"
distRows <- function(query, matrix, dist_func) {
  distances <- map(rownames(matrix), \(x) dist_func(matrix[x, ], query)) %>%
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
    result <- distRows(query[q, ], points, dist_func)
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


#' Plot dimensionality reduction results
#'
#' @description
#' Generic function for plotting dimensionality
#' reduction results, highlighting on specific column
#' Generates both 3d and 2d plots
plotDr <- function(to_plot, color_col, path, technique, labels, twod) {
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  if (technique == "pca") {
    prefix <- "PC"
  } else {
    prefix <- "V"
  }
  caption_str <- ifelse(labels$caption == "", "",
    glue("-{gsub(' ', '_', labels$caption)}")
  )
  if (twod) {
    biplot <- biplotCustom(to_plot,
      x = glue("{prefix}1"), y = glue("{prefix}2"),
      colour_column = color_col,
      labels = labels
    )
    mySaveFig(biplot, glue("{path}/{technique}_biplot-{color_col}{caption_str}.png"))
  } else {
    v1 <- to_plot[[glue("{prefix}1")]]
    v2 <- to_plot[[glue("{prefix}2")]]
    v3 <- to_plot[[glue("{prefix}3")]]
    plotly::plot_ly(to_plot,
      x = v1,
      y = v2,
      z = v3, color = ~ base::get(color_col),
      type = "scatter3d",
      marker = list(size = 5)
    ) %>% mySaveFig(., glue("{path}/{technique}_3d-{color_col}{caption_str}.html"))
  }
  return(biplot)
}

labelGen <- function(analysis_name, sample_name, caption) {
  return(list(
    title = glue("{analysis_name} of {sample_name}"),
    caption = caption
  ))
}
