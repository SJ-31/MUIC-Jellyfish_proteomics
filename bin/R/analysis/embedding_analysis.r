library("optparse")
library("glue")
parser <- OptionParser()
parser <- add_option(parser, c("-f", "--figure_path"))
parser <- add_option(parser, c("-s", "--sample_name"))
parser <- add_option(parser, c("-u", "--uniprot_data_dir"))
parser <- add_option(parser, "--technique")
parser <- add_option(parser, "--uniprot_embeddings")
parser <- add_option(parser, c("-m", "--mode")) # Dimensionality reduction technique
parser <- add_option(parser, "--protein_embedding_mode") # One of mean or sum
parser <- add_option(parser, c("-o", "--ontologizer_path"))
parser <- add_option(parser, c("-e", "--embeddings_path"))
parser <- add_option(parser, c("-r", "--r_source"))
parser <- add_option(parser, c("-p", "--python_source"))
parser <- add_option(parser, c("-c", "--combined_results"))
args <- parse_args(parser)
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))

# Gets all the necessary data
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))

switch(args$technique,
       "PCA" = {
         protein_path <- glue("{args$figure_path}/pca/protein_embeddings")
         go_path <- glue("{args$figure_path}/pca/go_embeddings")
         dr_func <- pcaAndJoin
       },
       "COSINE" = {
         protein_path <- glue("{args$figure_path}/cosine/protein_embeddings")
         go_path <- glue("{args$figure_path}/cosine/go_embeddings")
         dr_func <- cosinePcoaAndJoin
       },
       "UMAP" = {
         protein_path <- glue("{args$figure_path}/umap/protein_embeddings")
         go_path <- glue("{args$figure_path}/umap/go_embeddings")
         dr_func <- umapAndJoin
       },
       "TSNE" = {
         protein_path <- glue("{args$figure_path}/tsne/protein_embeddings")
         go_path <- glue("{args$figure_path}/tsne/go_embeddings")
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
