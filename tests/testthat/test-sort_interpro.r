library(testthat)
library(tidyverse)
source("./bin/sort_interpro.r")
pth <- "./tests/results/Unmatched"
real <- "./results/jellyfish//1-First_pass/Unmatched/"
args <- list(
  interpro_results = glue("{real}/InterPro/jellyfish_interpro.tsv"),
  eggnog_unmatched = glue("{real}/eggNOG/jellyfish_eggnog_unmatched.tsv"),
  output = glue("{pth}/Interpro/jellyfish_interpro_matched.tsv")
)
final_unmatched <- "testthat/output/interpro_unmatched"

results <- main(args)
output <- results$all
expect_equal(any(is.na(output$ID_method)), FALSE)
expect_equal(output$ProteinId, unique(output$ProteinId))
