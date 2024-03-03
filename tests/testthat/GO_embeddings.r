SAVE <- TRUE
MODE <- ""
args <- list(
  r_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/R",
  python_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/"
)
library(clv)
library(clValid)
source("./bin/R/GO_helpers.r")
source("./bin/R/DR_helpers.r")

fig_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/1-First_pass/Figures"
test_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/1-First_pass/C_indra_all.tsv"
embeddings_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/go_embedded.npz"
onto_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/ontologizer/"
sample_name <- "C_indra"
d <- goDataGlobal(
  uniprot_data_dir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa",
  onto_path = onto_path,
  sample_data = test_path,
  sample_name = sample_name
)


# Working with go embeddings
# Two sets of embeddings are being used here: "sample", which contains
# GO terms/proteins found from the proteins identified in the shotgun
# pipeline (i.e. from the C. indra data) and "all", which also
# contain GO terms/proteins of other taxa used for comparison

# Obtain embeddings from python
py$wanted_gos <- d$go_vec$all
py_run_string("wanted_gos = set(wanted_gos)")
all_embd_go <- py$loadEmbeddings(
  embeddings_path, "embds",
  py$wanted_gos
) %>%
  as_tibble() %>%
  t() %>%
  m2Tb(., first_col = "GO_IDs")
d$go_tb$all <- d$go_tb$all %>% dplyr::filter(GO_IDs %in% all_embd_go$GO_IDs)

if (!SAVE) {
  all_embd_prot_sum <- goEmbedding2Prot(d$prot_go_map$all, all_embd_go, sum) %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
  all_embd_prot_mean <- goEmbedding2Prot(d$prot_go_map$all, all_embd_go, mean) %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
} else {
  all_embd_prot_sum <- read_tsv("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/Analysis/all_embeddings_prot_sum.tsv") %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
  all_embd_prot_mean <- read_tsv("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/Analysis/all_embeddings_prot_mean.tsv") %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
}
sample_embd_prot_sum <- all_embd_prot_sum %>%
  dplyr::filter(ProteinId %in% d$sample_tb$ProteinId)
sample_embd_prot_mean <- all_embd_prot_mean %>%
  dplyr::filter(ProteinId %in% d$sample_tb$ProteinId)
sample_embd_go <- all_embd_go %>% filter(GO_IDs %in% d$go_vec$sample)
sample_go <- list(
  data = list(go = sample_embd_go),
  tb = dplyr::filter(
    d$go_tb$sample,
    is.element(
      GO_IDs,
      sample_embd_go$GO_IDs
    )
  ),
  color = c("sig_downloaded_db", "sig_id_w_open")
)
sample_protein <- list(
  data = list(
    sum = sample_embd_prot_sum,
    mean = sample_embd_prot_mean
  ),
  tb = d$sample_tb,
  color = c("Anno_method", "ID_method", "category")
)
all_go <- list(
  data = list(go = all_embd_go),
  tb = dplyr::filter(
    d$go_tb$all,
    is.element(
      GO_IDs,
      all_embd_go$GO_IDs
    )
  ),
  color = "taxon"
)
all_protein <- list(
  data = list(
    sum = all_embd_prot_sum,
    mean = all_embd_prot_mean
  ),
  tb = d$protein$all,
  color = "taxon"
)

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


# nearestInDim(toxin_embd$mean, 50, sample_embd_prot_mean,
#              euclideanDistance, "ProteinId")

#' Naming convention for saving figures
#' datatype_dataorigin_plottype-markerinfo
#' Figure path format is root/dr_type/data_type(GO|protein)

switch(MODE,
       "PCA" = {
         protein_path <- glue("{fig_path}/pca/protein_embeddings")
         go_path <- glue("{fig_path}/pca/go_embeddings")
         dr_func <- pcaAndJoin
       },
       "COSINE" = {
         protein_path <- glue("{fig_path}/cosine/protein_embeddings")
         go_path <- glue("{fig_path}/cosine/go_embeddings")
         dr_func <- cosinePcoaAndJoin
       },
       "UMAP" = {
         protein_path <- glue("{fig_path}/umap/protein_embeddings")
         go_path <- glue("{fig_path}/umap/go_embeddings")
         dr_func <- umapAndJoin
       },
       "TSNE" = {
         protein_path <- glue("{fig_path}/tsne/protein_embeddings")
         go_path <- glue("{fig_path}/tsne/go_embeddings")
         dr_func <- tsneSkAndJoin
       }
)
# toxin_pca_mean <- sample_pca$dr$mean$prcomp$x[names(toxin_proteins),] %>%
#   nearestInDim(query = ., data = sample_pca$dr$mean$prcomp, k = 50, dist_func = euclideanDistance)
# toxin_pca_sum <- sample_pca$dr$sum$prcomp$x[names(toxin_proteins),] %>%
#   nearestInDim(query = ., data = sample_pca$dr$sum$prcomp, k = 50, dist_func = euclideanDistance)

sample_go <- completeDR(sample_go, go_path,
                        "GO_IDs",
                        prefix = "sample", dr_func,
                        title = MODE
)
all_go <- completeDR(all_go, go_path,
                     "GO_IDs",
                     prefix = "all", dr_func,
                     title = MODE
)
sample_protein <- completeDR(sample_protein, protein_path,
                             "ProteinId",
                             prefix = "sample", dr_func,
                             title = MODE
)
all_protein <- completeDR(all_protein, protein_path,
                          "ProteinId",
                          prefix = "all", dr_func,
                          title = MODE
)


#
# Try to cluster embeddings
# Could you use the visualization to determine the number of groups?
# TODO:
dist_t <- sample_embd_prot_sum %>%
  t2Df(., "ProteinId") %>%
  dist()
clusters <- hclust(dist_t, method = "average")
plot(clusters, labels = FALSE)
cut <- cutree(clusters, h = 50)
# Look into dunn index to validate clusters
dunn(distance = dist_t, clusters = cut)
sil <- silhouette(dist = dist_t, x = cut)

scatt <- clv.Scatt(data = t2Df(sample_embd_prot_sum, "ProteinId"), clust = cut)
dens <- clv.DensBw(data = t2Df(sample_embd_prot_sum, "ProteinId"), clust = cut, scatt)
# Why is density infinity?
# The stack overflow post suggests it won't work with noise
# https://stackoverflow.com/questions/13177720/what-is-the-meaning-of-inf-in-s-dbw-output-in-r-commander

## clv.SDbw(scatt = scatt, dens = dens)
