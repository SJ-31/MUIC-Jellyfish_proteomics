library(tidyverse)
args <- list(r_source = "./bin/R")
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))

PATH <- "./results/C_indra_A"
run <- get_run("C_indra", TRUE, PATH)
alignments <- get_alignment_data(PATH)

combined <- run$first
mapping <- combined %>%
  filter(!is.na(MatchedPeptideIds)) %>%
  select(c(MatchedPeptideIds, ProteinId)) %>%
  separate_longer_delim( "MatchedPeptideIds", ";")

write_tsv(mapping, "peptides2proteinId.tsv")
