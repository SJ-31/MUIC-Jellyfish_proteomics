library("optparse")
library("glue")
parser <- OptionParser()
parser <- add_option(parser, c("-f", "--figure_path"))
parser <- add_option(parser, c("-s", "--sample_name"))
parser <- add_option(parser, c("-u", "--uniprot_data"))
parser <- add_option(parser, "--summed_embeddings")
parser <- add_option(parser, "--mean_embeddings")
parser <- add_option(parser, c("-m", "--mode"))
parser <- add_option(parser, c("-o", "--ontologizer_path"))
parser <- add_option(parser, c("-e", "--embeddings_path"))
parser <- add_option(parser, "--optimize",
  action = "store_true",
  default = FALSE
)
parser <- add_option(parser, c("-r", "--r_source"))
parser <- add_option(parser, c("-p", "--python_source"))
parser <- add_option(parser, c("-c", "--combined_results"))
args <- parse_args(parser)
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
fig_path <- args$figure_path

d <- goDataGlobal(
  uniprot_data_dir = args$uniprot_data,
  onto_path = args$ontologizer_path,
  sample_data = args$combined_results,
  sample_name = args$sample_name
)
py$wanted_gos <- d$go_vec$all
py_run_string("wanted_gos = set(wanted_gos)")
all_embd_go <- py$loadEmbeddings(
  args$embeddings_path, "embds",
  py$wanted_gos
) %>%
  as_tibble() %>%
  t() %>%
  m2Tb(., first_col = "GO_IDs")
d$go_tb$all <- d$go_tb$all %>% dplyr::filter(GO_IDs %in% all_embd_go$GO_IDs)
if (file.exists(args$summed_embeddings)) {
  all_embd_prot_sum <- read_tsv(args$summed_embeddings) %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
  all_embd_prot_mean <- read_tsv(args$mean_embeddings) %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
} else {
  all_embd_prot_sum <- goEmbedding2Prot(d$prot_go_map$all, all_embd_go, sum) %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
  write_tsv(all_embd_prot_sum, args$summed_embeddings)
  all_embd_prot_mean <- goEmbedding2Prot(d$prot_go_map$all, all_embd_go, mean) %>% dplyr::filter(ProteinId %in% d$protein$all$ProteinId)
  write_tsv(all_embd_prot_mean, args$mean_embeddings)
}

## Prepare embeddings
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

if (!args$optimize) {
  switch(args$mode,
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
      dr_func <- tsneAndJoin
    }
  )
  sample_go <- completeDR(sample_go, go_path,
    "GO_IDs",
    prefix = "sample", dr_func
  )
  all_go <- completeDR(all_go, go_path,
    "GO_IDs",
    prefix = "all", dr_func
  )
  sample_protein <- completeDR(sample_protein, protein_path,
    "ProteinId",
    prefix = "sample", dr_func
  )
  all_protein <- completeDR(all_protein, protein_path,
    "ProteinId",
    prefix = "all", dr_func
  )
} else {
  if (args$mode == "UMAP") {
    for (metric in c("cosine", "euclidean")) {
      for (n_neighbors in seq(from = 10, to = 15)) {
        for (min_dist in seq(from = 0.1, to = 0.5, by = 0.1)) {
          params <- list(umap = umap.defaults)
          params$umap <- umap.defaults
          params$umap$n_neighbors <- n_neighbors
          params$umap$min_dist <- min_dist
          params$umap$metric <- metric
          params$umap$n_components <- 3
          print(glue("Running umap with {n_neighbors} neighbors, min_dist: {min_dist}, {metric} distance"))
          path <- glue("{fig_path}/umap_optimization/{n_neighbors}nn_{min_dist}mind_{metric}")
          label <- labelGen("UMAP", "test", glue("n neighbors: {n_neighbors}, min dist: {min_dist}, metric: {metric}"))
          completeDR(
            dr_data = sample_protein,
            fig_dir = path,
            prefix = "sample",
            join_on = "ProteinId",
            dR = umapAndJoin,
            params = params,
            label = label
          )
        }
      }
    }
  } else if (args$mode == "TSNE") {
    for (metric in c("cosine", "euclidean")) {
      for (p in seq(from = 5, to = 50, by = 5)) {
        params <- list(
          perplexity = p,
          metric = metric
        )
        print(glue("Running tsne with perplexity {p}, {metric} distance"))
        path <- glue("{fig_path}/tsne_optimization/{p}p_{metric}")
        label <- labelGen("T-SNE", "test", glue("perplexity: {p}, metric: {metric}"))
        completeDR(
          dr_data = sample_protein,
          fig_dir = path,
          prefix = "sample",
          join_on = "ProteinId",
          dR = tsneSkAndJoin,
          params = params,
          label = label
        )
      }
    }
  }
}
