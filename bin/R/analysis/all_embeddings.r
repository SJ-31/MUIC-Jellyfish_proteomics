library("optparse")
library("glue")
parser <- OptionParser()
parser <- add_option(parser, c("-f", "--figure_path"))
parser <- add_option(parser, c("-s", "--sample_name"))
parser <- add_option(parser, c("-u", "--uniprot_data"))
parser <- add_option(parser, c("-m", "--mode"))
parser <- add_option(parser, "--protein_embd_mode")
parser <- add_option(parser, c("-o", "--ontologizer_path"))
parser <- add_option(parser, c("-e", "--embeddings_path"))
parser <- add_option(parser, c("-e", "--embeddings_path"))
parser <- add_option(parser, "--optimize",
                     action = "store_true",
                     default = FALSE
)
parser <- add_option(parser, c("-r", "--r_source"))
parser <- add_option(parser, c("-c", "--combined_results"))
args <- parse_args(parser)
source(glue("{args$r_source}/GO.r"))
source(glue("{args$r_source}/DR.r"))
fig_path <- args$figure_path

source("./bin/R/analysis/prepare_embeddings.r")

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
  sample_go <- completeDR(d$go_vec$sample, go_path,
                          "GO_IDs",
                          prefix = "sample", dr_func
  )
  all_go <- completeDR(d$go_vec$all, go_path,
                       "GO_IDs",
                       prefix = "all", dr_func
  )
  sample_protein <- completeDR(d$protein$sample, protein_path,
                               "ProteinId",
                               prefix = "sample", dr_func
  )
  all_protein <- completeDR(d$protein$all, protein_path,
                            "ProteinId",
                            prefix = "all", dr_func
  )
}
