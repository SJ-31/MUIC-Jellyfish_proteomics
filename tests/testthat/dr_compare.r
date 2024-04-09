library(optparse)
library(ggpubr)
library(glue)
if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-f", "--figure_path"))
  parser <- add_option(parser, c("-s", "--sample_name"))
  parser <- add_option(parser, c("-u", "--uniprot_data_dir"))
  parser <- add_option(parser, "--technique")
  parser <- add_option(parser, c("-o", "--ontologizer_path"))
  parser <- add_option(parser, "--metric")
  parser <- add_option(parser, "--uniprot_embeddings")
  parser <- add_option(parser, "--protein_embedding_mode") # One of mean or sum
  parser <- add_option(parser, c("-e", "--embeddings_path"))
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-p", "--python_source"))
  parser <- add_option(parser, c("-c", "--combined_results"))
  args <- parse_args(parser)
}
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
metric <- args$metric # One of euclidean, cosine

sample_protein$color <- "category"
results <- list()
min_dist <- 0.1
if (args$technique == "UMAP") {
  for (n_neighbors in seq(from = 5, to = 50, by = 10)) {
    params <- list(umap = umap.defaults)
    params$umap <- umap.defaults
    params$umap$n_neighbors <- n_neighbors
    params$umap$min_dist <- min_dist
    params$umap$metric <- metric
    params$umap$n_components <- 3
    print(glue("Running umap with {n_neighbors} neighbors, min_dist: {min_dist}, {metric} distance"))
    path <- glue("{args$figure_path}/umap_optimization/{n_neighbors}nn_{min_dist}mind_{metric}")
    label <- labelGen("UMAP", "test", glue("n neighbors: {n_neighbors}, min dist: {min_dist}, metric: {metric}"))
    dr <- drWrapper(
      dr_data = sample_protein,
      outdir = path,
      name = "sample",
      join_on = "ProteinId",
      technique = "UMAP",
      params = params
    )
    bp <- plotDr(dr$to_plot, sample_protein$color, "UMAP", label)
    results <- append(results, bp)
  }
} else if (args$technique == "TSNE") {
  for (p in seq(from = 5, to = 50, by = 5)) {
    params <- list(
      perplexity = p,
      metric = metric
    )
    print(glue("Running tsne with perplexity {p}, {metric} distance"))
    path <- glue("{args$figure_path}/tsne_optimization/{p}p_{metric}")
    label <- labelGen("T-SNE", "test", glue("perplexity: {p}, metric: {metric}"))
    dr <- drWrapper(
      dr_data = sample_protein,
      outdir = path,
      name = "sample",
      join_on = "ProteinId",
      technique = "TSNE",
      params = params
    )
    bp <- plotDr(dr$to_plot, sample_protein$color, "TSNE", label)
    results <- append(results, bp)
  }
}
comparison <- ggarrange(plotlist = results, common.legend = TRUE)
ggsave(glue("{args$figure_path}/{args$technique}_{metric}_comparison.png"),
       comparison,
       width = 30, height = 25
)
