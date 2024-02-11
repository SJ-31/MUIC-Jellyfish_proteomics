# Helper functions for analyzing GO embeddings and performing DR

args <- list(r_source = "./bin")
library(tidyverse)
library(glue)
library(umap)
library(tsne)
library(reticulate)
source_python(glue::glue("{args$r_source}/a2v.py"))
source(glue::glue("{args$r_source}/rrvgo_modified.r"))


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

biplotCustom <- function(ordination_tb, colour_column, x, y, palette) {
  # Convenience function for plotting ordination results
  if (missing(x) || missing(y)) {
    y <- "PC1"
    x <- "PC2"
  }
  if (missing(palette)) {
    p <- "Redmonder::qPBI"
  } else {
    p <- palette
  }
  plotted <- ggplot(ordination_tb, aes(
    x = .data[[x]], y = .data[[y]],
    colour = .data[[colour_column]]
  )) +
    geom_point(size = 2.5) +
    scale_color_paletteer_d(p) +
    theme_bw()
  return(plotted)
}

goEmbedding2Prot <- function(protein_map, embedding_tb, combine_func) {
  # Combine GO terms assigned to each protein in protein_map using
  # combine_func
  purrr::map(names(protein_map), ~{
    dplyr::filter(embedding_tb, GO_IDs %in% protein_map[[.x]]) %>%
      reframe(across(is.numeric, combine_func)) %>%
      mutate(ProteinId = .x, .before = V1)
  }) %>% bind_rows()
}

# Wrapper for tsne on "data", removing the non-numeric join_col
tsneAndJoin <- function(data, join_col) {
  data_only <- data %>% dplyr::select(-!!join_col)
  tsne_tb <- tsne(data_only, k = 3)
  return(tsne_tb)
}

# Perform umap on "data",
umapAndJoin <- function(data, join_col) {
  umap_params <- umap.defaults
  umap_params$n_neighbors <- 10
  data_only <- data %>% dplyr::select(-!!join_col)
  umap_params$n_components <- 3
  umap <- umap(data_only, config = umap_params)
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
plotDr <- function(dr_result, join_tb, join_col, color_col, path, name) {
  prefix <- "V"
  if (pluck_exists(dr_result, "prcomp")) {
    prefix <- "PC"
    dr_result <- dr_result$prcomp[["x"]]
  } else if (class(dr_result) == "wcmdscale") {
    dr_result <- dr_result$points
    colnames(dr_result) <- colnames(dr_result) %>%
      map_chr(., \(x) gsub("Dim", "V", x))
  }
  plot <- dr_result %>%
    as_tibble() %>%
    mutate(!!join_col := join_tb[[join_col]]) %>%
    inner_join(., join_tb, by = join_by(!!join_col))
  biplotCustom(plot,
               x = glue("{prefix}1"), y = glue("{prefix}2"),
               colour_column = color_col
  ) %>%
    mySaveFig(., glue("{path}/{name}_biplot-{color_col}.png"))
  plot_ly(plot,
          x = ~get(glue("{prefix}1")),
          y = ~get(glue("{prefix}2")),
          z = ~get(glue("{prefix}3")), color = ~get(color_col)
  ) %>%
    mySaveFig(., glue("{path}/{name}_3d-{color_col}.html"))
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
completeDR <- function(dr_data, fig_dir, join_on, prefix, dR) {
  for (n in names(dr_data$data)) {
    data <- purrr::pluck(dr_data, "data", n)
    dr_data$dr[[n]] <- dR(data, join_on)
    purrr::map(dr_data$color, \(x) plotDr(dr_data$dr[[n]],
                                          color_col = x,
                                          join_tb = dr_data$tb,
                                          join_col = join_on,
                                          path = paste0(fig_dir, "_", n),
                                          name = paste0(prefix, "_", x)))
  }
  return(dr_data)
}
