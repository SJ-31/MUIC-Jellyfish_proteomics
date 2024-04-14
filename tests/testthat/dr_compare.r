library(optparse)
library(ggpubr)
library(glue)
if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-f", "--figure_path"))
  parser <- add_option(parser, c("-c", "--combined_results"))
  parser <- add_option(parser, c("-s", "--sample_name"))
  parser <- add_option(parser, c("-d", "--dist_path"))
  parser <- add_option(parser, "--technique")
  parser <- add_option(parser, "--metric")
  parser <- add_option(parser, c("-e", "--embedding_path"))
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-p", "--python_source"))
  args <- parse_args(parser)
}
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
metric <- args$metric # One of euclidean, cosine


e <- embeddingData(
  args$combined_results,
  args$sample_name,
  args$embedding_path,
  args$dist_path,
)
ARRANGE <- TRUE

e$color <- "category"
results <- list()
min_dist <- 0.1

if (args$technique == "pca" || args$technique == "pcoa") {
  path <- glue("{args$figure_path}/{args$technique}")
  dr <- drWrapper(
    dr_data = e,
    outdir = path,
    name = "sample",
    join_on = "ProteinId",
    technique = args$technique
  )
  label <- labelGen(args$technique, "test", "")
  bp <- plotDr(dr$to_plot, e$color, path, args$technique, label, TRUE)
  ARRANGE <- FALSE
}


if (args$technique == "umap") {
  for (n_neighbors in seq(from = 5, to = 50, by = 10)) {
    PARAMS$umap$n_neighbors <- n_neighbors
    PARAMS$umap$min_dist <- min_dist
    PARAMS$metric <- metric
    PARAMS$umap$n_components <- 2
    print(glue("Running umap with {n_neighbors} neighbors, min_dist: {min_dist}, {metric} distance"))
    path <- glue("{args$figure_path}/umap_optimization/{n_neighbors}nn_{min_dist}mind_{metric}")
    label <- labelGen("UMAP", "test", glue("n neighbors: {n_neighbors}, min dist: {min_dist}, metric: {metric}"))
    dr <- drWrapper(
      dr_data = e,
      outdir = path,
      name = "sample",
      join_on = "ProteinId",
      technique = "umap"
    )
    bp <- plotDr(dr$to_plot, e$color, path, "umap", label, TRUE)
    results <- append(results, bp)
  }
} else if (args$technique == "tsne") {
  for (p in seq(from = 5, to = 50, by = 5)) {
    PARAMS$metric <- metric
    PARAMS$tsne$perplexity <- p
    print(glue("Running tsne with perplexity {p}, {metric} distance"))

    path <- glue("{args$figure_path}/tsne_optimization/{p}p_{metric}")
    label <- labelGen("T-SNE", "test", glue("perplexity: {p}, metric: {metric}"))
    dr <- drWrapper(
      dr_data = e,
      outdir = path,
      name = "sample",
      join_on = "ProteinId",
      technique = "tsne"
    )
    bp <- plotDr(dr$to_plot, e$color, path, "tsne", label, TRUE)
    results <- append(results, bp)
  }
}


if (ARRANGE) {
  comparison <- ggarrange(plotlist = results)
  ggsave(glue("{args$figure_path}/{args$technique}_{metric}_comparison.png"),
         comparison,
         width = 30, height = 25
  )
}
