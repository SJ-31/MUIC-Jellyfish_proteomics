SAVE <- TRUE
MODE <- ""
args <- list(
  r_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/R",
  python_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/",
  uniprot_embeddings = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/uniprot_embeddings_mean.tsv",
  uniprot_data_dir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa",
  figure_path = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/testthat/output/figures",
  combined_results = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/1-First_pass/C_indra_all.tsv",
  embeddings_path = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/go_embedded.npz",
  ontologizer_path = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/ontologizer/",
  sample_name <- "C_indra"
)
library(clv)
library(ggpubr)
library(clValid)
source("./bin/R/GO_helpers.r")
source("./bin/R/DR_helpers.r")

protein_embd_mode <- "mean" # One of mean or sum

source("../../bin/R/analysis/prepare_embeddings.r")


# A list of k nearest proteins to the query proteins, using different distance metrics
nearest <- list()
toxin_proteins <- getToxinProteins(d$prot_go_map$sample)
toxin_embd <- list(
  mean = dplyr::filter(
    sample_embd_prot_mean,
    ProteinId %in% names(toxin_proteins)
  ),
  sum =
    dplyr::filter(
      sample_embd_prot_sum,
      ProteinId %in% names(toxin_proteins)
    )
) %>% lapply(., \(x) column_to_rownames(x, var = "ProteinId") %>% as.matrix())


#' Naming convention for saving figures
#' datatype_dataorigin_plottype-markerinfo
#' Figure path format is root/dr_type/data_type(GO|protein)

# toxin_pca_mean <- sample_pca$dr$mean$prcomp$x[names(toxin_proteins),] %>%
#   nearestInDim(query = ., data = sample_pca$dr$mean$prcomp, k = 50, dist_func = euclideanDistance)
# toxin_pca_sum <- sample_pca$dr$sum$prcomp$x[names(toxin_proteins),] %>%
#   nearestInDim(query = ., data = sample_pca$dr$sum$prcomp, k = 50, dist_func = euclideanDistance)

#
# Try to cluster embeddings
# Could you use the visualization to determine the number of groups?
# TODO:
## dist_t <- sample_embd_prot_sum %>%
##   t2Df(., "ProteinId") %>%
##   dist()
## clusters <- hclust(dist_t, method = "average")
## plot(clusters, labels = FALSE)
## cut <- cutree(clusters, h = 50)
## # Look into dunn index to validate clusters
## dunn(distance = dist_t, clusters = cut)
## sil <- silhouette(dist = dist_t, x = cut)

## scatt <- clv.Scatt(data = t2Df(sample_embd_prot_sum, "ProteinId"), clust = cut)
## dens <- clv.DensBw(data = t2Df(sample_embd_prot_sum, "ProteinId"), clust = cut, scatt)
# Why is density infinity?
# The stack overflow post suggests it won't work with noise
# https://stackoverflow.com/questions/13177720/what-is-the-meaning-of-inf-in-s-dbw-output-in-r-commander

## clv.SDbw(scatt = scatt, dens = dens)

metric <- "euclidean"
min_dist <- 0.1
results <- list()
sample_protein$color <- "category"
for (n_neighbors in seq(from = 10, to = 15)) {
  ## for (min_dist in seq(from = 0.1, to = 0.5, by = 0.1)) {
  params <- list(umap = umap.defaults)
  params$umap <- umap.defaults
  params$umap$n_neighbors <- n_neighbors
  params$umap$min_dist <- min_dist
  params$umap$metric <- metric
  params$umap$n_components <- 3
  print(glue("Running umap with {n_neighbors} neighbors, min_dist: {min_dist}, {metric} distance"))
  path <- glue("{args$figure_path}/umap_optimization/{n_neighbors}nn_{min_dist}mind_{metric}")
  label <- labelGen("UMAP", "test", glue("n neighbors: {n_neighbors}, min dist: {min_dist}, metric: {metric}"))
  dr <- completeDR(
    dr_data = sample_protein,
    fig_dir = path,
    prefix = "sample",
    join_on = "ProteinId",
    dR = umapAndJoin,
    params = params,
    label = label
  )
  results <- append(results, dr$biplots)
  ## }
}

comparison <- ggarrange(plotlist = results, common.legend = TRUE)
ggsave("~/test_compare.png", comparison, width = 20, height = 15)
