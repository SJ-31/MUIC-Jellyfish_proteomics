library("optparse")
library("glue")
parser <- OptionParser()
parser <- add_option(parser, c("-f", "--figure_path"))
parser <- add_option(parser, c("-s", "--sample_name"))
parser <- add_option(parser, c("-u", "--uniprot_data"))
parser <- add_option(parser, c("-t", "--technique"))
parser <- add_option(parser, c("-e", "--embeddings_path"))
parser <- add_option(parser, c("-r", "--r_source"))
parser <- add_option(parser, c("-r", "--python_source"))
parser <- add_option(parser, "--results_prefix")
parser <- add_option(parser, "--compare", action = "store_true")
parser <- add_option(parser, c("-c", "--combined_results"))
args <- parse_args(parser)
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))

source("./bin/R/analysis/prepare_embeddings.r")
if (args$compare) {
  e <- embeddingData(args$combined_results,
                     args$sample_name,
                     args$embedding_path,
                     args$dist_path, FALSE,
                     comparison_meta = args$uniprot_data)
} else {
  e <- embeddingData(args$combined_results,
                     args$sample_name,
                     args$embedding_path,
                     args$dist_path, TRUE)
}


prefix <- ifelse(!args$results_prefix, args$sample_name,
                 args$results_prefix)
result <- drWrapper(e, "ProteinId",
                    args$figure_path,
                    prefix,
                    args$technique)
title_str <- ifelse(args$compare,
                    glue("{args$sample_name} comparison"),
                    glue("{args$sample_name} sample"))

for (color in data$color) {
  label <- labelGen(args$technique, title_str)
  plotDr(result$to_plot, args$figure_path, args$technique,
         labels = label, twod_only = FALSE)
}
