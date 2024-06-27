# Script for calculating GO embeddings of database proteins
library(glue)
if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-r", "--r_source"), type = "character")
  parser <- add_option(parser, c("-p", "--python_source"), type = "character")
  parser <- add_option(parser, c("-u", "--uniprot_data"), type = "character")
  parser <- add_option(parser, c("-e", "--embeddings_path"), type = "character")
  parser <- add_option(parser, c("-o", "--output"), type = "character")
  args <- parse_args(parser)
}
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
u <- get_uniprot(args$uniprot_data)

# Obtain embeddings from python
py$wanted_gos <- u$go_vec
py_run_string("wanted_gos = set(wanted_gos)")
all_embd_go <- py$loadEmbeddings(
  args$embeddings_path, "embds",
  py$wanted_gos
) %>%
  as_tibble() %>%
  t() %>%
  m2tb(., first_col = "GO_IDs")

all_embd_prot_sum <- goEmbedding2Prot(u$map, all_embd_go, sum)
all_embd_prot_mean <- goEmbedding2Prot(u$map, all_embd_go, mean)
write_tsv(all_embd_prot_mean, glue("{args$output}/uniprot_embeddings_mean.tsv"))
write_tsv(all_embd_prot_sum, glue("{args$output}/uniprot_embeddings_sum.tsv"))
