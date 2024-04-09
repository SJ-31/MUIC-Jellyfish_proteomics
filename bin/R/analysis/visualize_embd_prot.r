library("optparse")
library("glue")
parser <- OptionParser()
parser <- add_option(parser, c("-f", "--figure_path"))
parser <- add_option(parser, c("-s", "--sample_name"))
parser <- add_option(parser, c("-u", "--uniprot_data"))
parser <- add_option(parser, c("-t", "--technique"))
parser <- add_option(parser, c("-o", "--ontologizer_path"))
parser <- add_option(parser, c("-e", "--embeddings_path"))
parser <- add_option(parser, c("-r", "--r_source"))
parser <- add_option(parser, c("-r", "--python_source"))
parser <- add_option(parser, c("-c", "--combined_results"))
args <- parse_args(parser)
source(glue("{args$r_source}/GO.r"))
source(glue("{args$r_source}/DR.r"))
fig_path <- args$figure_path
args$embd_type <- "protein"

source("./bin/R/analysis/prepare_embeddings.r")
