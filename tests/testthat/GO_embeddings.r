SAVE <- TRUE
PCA <- TRUE
UMAP <- FALSE
TSNE <- FALSE
COSINE <- FALSE
library(plotly)
source("./bin/GO.r")
source("./bin/DR.r")

fig_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/1-First_pass/Figures"
test_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/1-First_pass/jellyfish_all.tsv"
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
all_embd_go <- py$loadEmbeddings(embeddings_path, "embds",
                                 py$wanted_gos) %>%
  as_tibble() %>%
  t() %>%
  m2Tb(., first_col = "GO_IDs")
d$go_tb$all <- d$go_tb$all %>% dplyr::filter(GO_IDs %in% all_embd_go$GO_IDs)

if (!SAVE) {
  all_embd_prot_sum <- goEmbedding2Prot(d$prot_go_map$all, all_embd_go, sum)
  all_embd_prot_mean <- goEmbedding2Prot(d$prot_go_map$all, all_embd_go, mean)
} else {
  all_embd_prot_sum <- read_tsv("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/Analysis/all_embeddings_prot_sum.tsv")
  all_embd_prot_mean <- read_tsv("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/Analysis/all_embeddings_prot_mean.tsv")
}
sample_embd_prot_sum <- all_embd_prot_sum %>%
  dplyr::filter(ProteinId %in% d$sample_tb$ProteinId)
sample_embd_prot_mean <- all_embd_prot_mean %>%
  dplyr::filter(ProteinId %in% d$sample_tb$ProteinId)
sample_embd_go <- all_embd_go %>% filter(GO_IDs %in% d$go_vec$sample)
sample_go <- list(data = list(go = sample_embd_go),
                  tb = dplyr::filter(d$go_tb$sample,
                                     is.element(GO_IDs,
                                                sample_embd_go$GO_IDs)),
                  color = c("sig_downloaded_db", "sig_id_w_open"))
sample_protein <- list(data = list(sum = sample_embd_prot_sum,
                                   mean = sample_embd_prot_mean),
                       tb = d$sample_tb,
                       color = c("Anno_method", "ID_method"))
all_go <- list(data = list(go = all_embd_go),
               tb = dplyr::filter(d$go_tb$all,
                                  is.element(GO_IDs,
                                             all_embd_go$GO_IDs)),
               color = "taxon")
all_protein <- list(data = list(sum = all_embd_prot_sum,
                                mean = all_embd_prot_mean),
                    tb = d$protein$all,
                    color = "taxon")

# Label known and unknown proteins for coloring, using header information
toxins <- markMatch(d$sample_tb, "header", "[tT]oxin", "toxin")
others <- markMatch(d$sample_tb, "header", "unknown", "unknown")

header_queries <- list(
  unknown = "unknown|uncharacterized|predicted|unnamed",
  venom_component = "toxin|porin",
  transport = "pump|transporter|transporting"
)

#' Naming convention for saving figures
#' datatype_dataorigin_plottype-markerinfo
#' Figure path format is root/dr_type/data_type(GO|protein)

#' PCA
#'
#' @description
#' All the pca computations and plotting
if (PCA) {
  pca_protein_path <- glue("{fig_path}/pca/protein_embeddings")
  pca_go_path <- glue("{fig_path}/pca/go_embeddings")
  sample_go <- completeDR(sample_go, pca_go_path,
                          "GO_IDs", prefix = "sample", pcaAndJoin)
  xvals <- sample_go$dr$go$prcomp$x

  reconstructed <- pcaReconstruct(sample_go$dr$go$prcomp, 2)
  reconstructed <- pcaReconstruct(sample_go$dr$go$prcomp)
  original <- sample_go$data$go %>%
    dplyr::select(-GO_IDs) %>%
    as.matrix() %>%
    sweep(., 2, colMeans(.))
  mean(as.matrix(dplyr::select(sample_go$data$go, -GO_IDs)) - reconstructed)
  assertthat::are_equal(xvals, original %*% sample_go$dr$go$prcomp$rotation)

  all_go <- completeDR(all_go, pca_go_path,
                       "GO_IDs", prefix = "all", pcaAndJoin)
  sample_pca <- completeDR(sample_protein, pca_protein_path,
                           "ProteinId",
                           prefix = "sample",
                           pcaAndJoin)
  all_pca <- completeDR(all_protein, pca_protein_path,
                        "ProteinId", prefix = "all", pcaAndJoin)
}


#' Cosine distance/dissimilarity
#'
#' @description
#' Computations and plotting using cosine dissimilarity
if (COSINE) {
  cosine_sample <- list()
  cosine_sample$pcoa <- sample_embd_prot_sum %>%
    cosineDissimilarity("ProteinId") %>%
    pcoaWithTb(distances = ., tb = d$sample_tb, join_on = "ProteinId")
  cosine_sample$biplot$anno <- cosine_sample$pcoa$points %>%
    biplotCustom(., colour_column = "Anno_method", x = "Dim1", y = "Dim2")

  cosine_all <- list()

}


#' UMAP
#'
#' @description
#' Computations and plotting for UMAP
if (UMAP) {
  umap_protein_path <- glue("{fig_path}/umap/protein_embeddings")
  umap_go_path <- glue("{fig_path}/umap/go_embeddings")
  sample_go_umap <- completeDR(sample_go, umap_go_path,
                               "GO_IDs",
                               prefix = "sample", umapAndJoin)
  all_go_umap <- completeDR(all_go, umap_go_path,
                            "GO_IDs", prefix = "all", umapAndJoin)
  sample_protein_umap <- completeDR(sample_protein, umap_protein_path,
                                    "ProteinId",
                                    prefix = "sample", umapAndJoin)
  all_protein_umap <- completeDR(all_protein, umap_protein_path,
                                 "ProteinId",
                                 prefix = "all", umapAndJoin)
}


#' t-SNE
#'
#' @description
#' Computations and plotting for t-SNE
if (TSNE) {
  tsne_protein_path <- glue("{fig_path}/tsne/protein_embeddings")
  tsne_go_path <- glue("{fig_path}/tsne/go_embeddings")
  sample_go_tsne <- completeDR(sample_go, tsne_go_path,
                               "GO_IDs",
                               prefix = "sample", tsneAndJoin)
  all_go_tsne <- completeDR(all_go, tsne_go_path,
                            "GO_IDs", prefix = "all", tsneAndJoin)
  sample_protein_tsne <- completeDR(sample_protein, tsne_protein_path,
                                    "ProteinId",
                                    prefix = "sample", tsneAndJoin)
  all_protein_tsne <- completeDR(all_protein, tsne_protein_path,
                                 "ProteinId",
                                 prefix = "all", tsneAndJoin)
}

#' Cosine pcoa
#'
#' @description
#' Computations and plotting for PcOA using cosine distance
if (COSINE) {

}


