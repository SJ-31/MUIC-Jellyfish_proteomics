library(testthat)
library(glue)
bin <- "./bin/R"
source(glue("{bin}/unify_groups.r"))
path <- "./results//C_indra_A/1-First_pass//Combined/intersected_searches.tsv"
path <- "~/Downloads/testdata.tsv"
m <- main("standard", "G", path)
## write_tsv(m, "./results/jellyfish/1-First_pass//Combined//unified_groups.tsv")
