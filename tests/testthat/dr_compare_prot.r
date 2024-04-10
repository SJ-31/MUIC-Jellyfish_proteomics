source("./bin/R/GO_helpers.r")
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  embd_type = "protein",
  sample_name = "C_indra",
  combined_results = glue("{wd}/results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv"),
  sample_embd = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/embeddings.hdf5"),
  sample_embd_dist = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"
  ),
  figure_path = glue("{wd}//tests/testthat/protein_dr_compare"))
source("./bin/R/DR_helpers.r")

source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
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
