library(testthat)
library(tidyverse)
source("./bin/sort_interpro.r")
pth <- "./tests/results/Unmatched"
args <- list(
  interpro_results = glue("{pth}/InterPro/nO_D_interpro.tsv"),
  eggnog_unmatched = glue("{pth}/eggNOG/nO_D_eggnog_unmatched.tsv")
)
anno_output <- "tests/nf-test-out/sort_interpro/interpro_anno-nO_D.tsv"
meta_output <- "tests/nf-test-out/sort_interpro/interpro_meta-nO_D.tsv"
final_unmatched <- "testthat/output/interpro_unmatched"

output <- main(args)
anno <- output$anno_df
meta <- output$meta_df
expect_equal(any(is.na(meta$ID_method)), FALSE)
