library(tidyverse)
library(testthat)
source("./bin/sort_eggnog.r")
pth <- "./tests/results/Unmatched"
args <- list(
  blast = glue("{pth}/eggNOG/nO_D_unmatched.tsv"),
  peptides = glue("{pth}/unmatched_peptides.tsv"),
  annotations = glue("{pth}/eggNOG/nO_D.emapper.annotations"),
  seeds = glue("{pth}/eggNOG/nO_D.emapper.seed_orthologs")
)
m <- main(args)
anno <- m$anno_df
meta <- m$meta_df
expect_equal(any(is.na(meta$ID_method)), FALSE)
u <- m$unmatched$tsv
expect_equal(any(is.na(u$ID_method)), FALSE)
expect_equal(meta$ProteinId, unique(meta$ProteinId))
