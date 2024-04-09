library("optparse")
library("glue")
parser <- OptionParser()
parser <- add_option(parser, c("-f", "--figure_path"))
parser <- add_option(parser, c("-s", "--sample_name"))
parser <- add_option(parser, c("-u", "--uniprot_data"))
parser <- add_option(parser, c("-t", "--technique"))
parser <- add_option(parser, "--protein_embd_mode")
parser <- add_option(parser, c("-o", "--ontologizer_path"))
parser <- add_option(parser, c("-r", "--r_source"))
parser <- add_option(parser, c("-r", "--python_source"))
parser <- add_option(parser, c("-c", "--combined_results"))
args <- parse_args(parser)
source(glue("{args$r_source}/GO.r"))
source(glue("{args$r_source}/DR.r"))
fig_path <- args$figure_path
args$embd_type <- "GO"

source("./bin/R/analysis/prepare_embeddings.r")

protein_path <- glue("{fig_path}/protein")
go_path <- glue("{fig_path}/go")

sample_go <- drWrapper(d$go_vec$sample, go_path,
                       "GO_IDs",
                       sample_name = glue("{args$sample_name}"),
                       args$technique
)
all_go <- drWrapper(d$go_vec$all, go_path,
                    "GO_IDs",
                    sample_name = glue("{args$sample_name} compare"),
                    args$technique
)
sample_protein <- drWrapper(d$protein$sample, protein_path,
                            "ProteinId",
                            sample_name = "sample",
                            args$technique
)
all_protein <- drWrapper(d$protein$all, protein_path,
                         "ProteinId",
                         sample_name = "all",
                         args$technique
)
