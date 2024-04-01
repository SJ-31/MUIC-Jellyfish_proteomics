args <- list(
  r_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/R",
  python_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/",
  uniprot_embeddings = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/uniprot_embeddings_mean.tsv",
  uniprot_data_dir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa",
  figure_path = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/testthat/output/figures",
  combined_results = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/1-First_pass/C_indra_all.tsv",
  embeddings_path = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/go_embedded.npz",
  ontologizer_path = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/ontologizer/",
  sample_name = "C_indra",
  protein_embd_mode = "mean" # One of mean or sum
)
library(clv)
library(fpc)
library(clValid)
library(leiden)
library(igraph)
source("./bin/R/GO_helpers.r")
source("./bin/R/DR_helpers.r")


SAMPLE_ONLY <- TRUE
source("./bin/R/analysis/prepare_embeddings.r")


# A list of k nearest proteins to the query proteins, using different distance metrics
# nearest <- list()
# toxin_proteins <- getToxinProteins(d$prot_go_map$sample)
# toxin_embd <- dplyr::filter(sample_prot_embd,
#                             ProteinId %in% names(toxin_proteins)) %>%
#   column_to_rownames(., var = "ProteinId") %>%
#   as.matrix()


#' Naming convention for saving figures
#' datatype_dataorigin_plottype-markerinfo
#' Figure path format is root/dr_type/data_type(GO|protein)

# toxin_pca_mean <- sample_pca$dr$mean$prcomp$x[names(toxin_proteins),] %>%
#   nearestInDim(query = ., data = sample_pca$dr$mean$prcomp, k = 50, dist_func = euclideanDistance)
# toxin_pca_sum <- sample_pca$dr$sum$prcomp$x[names(toxin_proteins),] %>%
#   nearestInDim(query = ., data = sample_pca$dr$sum$prcomp, k = 50, dist_func = euclideanDistance)


## and leiden

