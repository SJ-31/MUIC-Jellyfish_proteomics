library(optparse)
library(glue)
PLOT <- FALSE
if (PLOT) {
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
} else {
  args <- list(r_source = "./bin/R",
               python_source = "./bin")
}
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))

if (PLOT) {
  metric <- args$metric # One of euclidean, cosine

  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    args$embedding_path,
    args$dist_path,
  )

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
}

if (grepl("Bio_SDD", getwd())) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests"
} else {
  wd <- "/home/shannc/workflow/tests"
}
library(paletteer)
outdir <- glue("{wd}/testthat/output/protein_dr")
pca_pcoa <- list.files(outdir, full.names = TRUE) %>%
  keep(grepl("pca|pcoa", .)) %>%
  lapply(\(x) {
    list.files(x, recursive = TRUE, full.names = TRUE) %>% keep(grepl(".*_explained.tsv", .))
  }) %>%
  unlist()
max_v <- 30 # For better visualization, don't keep all the variables
tbs <- pca_pcoa %>% lapply(., \(x) {
  description <- str_split_1(x, "/")[11] %>% str_split_1("_")
  technique <- description[1]
  model <- description[3]
  tb <- read_tsv(x) %>%
    mutate(technique = technique,
           model = model) %>%
    slice(1:max_v)
  return(tb)
})
pca <- tbs %>%
  keep(\(x) x$technique[1] == "pca") %>%
  bind_rows()
pcoa <- tbs %>%
  keep(\(x) x$technique[1] == "pcoa") %>%
  bind_rows()

pca_plot <- pca %>% ggplot(aes(x = PC, y = Variance_explained, fill = model)) +
  scale_fill_paletteer_d("feathers::eastern_rosella") +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Variance explained (%)") +
  scale_x_continuous(breaks = seq(1, max_v, by = 1)) +
  geom_line(aes(x = PC, y = Cum_sum_variance_explained, color = model)) +
  scale_color_paletteer_d("feathers::eastern_rosella")
ggsave(glue("{outdir}/pca_explained.png"), pca_plot)

pcoa_plot <- pcoa %>% ggplot(aes(x = PCo, y = Percent_explained, fill = model)) +
  scale_fill_paletteer_d("fishualize::Lepomis_megalotis") +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Percent explained (%)") +
  scale_x_continuous(breaks = seq(1, max_v, by = 1)) +
  geom_line(aes(x = PCo, y = Cum_sum_percent_explained, color = model)) +
  scale_color_paletteer_d("fishualize::Lepomis_megalotis")

ggsave(glue("{outdir}/pcoa_explained.png"), pcoa_plot)


