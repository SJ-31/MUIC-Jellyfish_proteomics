library(testthat)
library(glue)
bin <- "./bin"
source(glue("{bin}/unify_groups.r"))
path <- "./tests/testthat//output//intersected_searches.tsv"
m <- main("standard", "G", path)
write_tsv(m, "./results/jellyfish//1-First_pass//Combined//unified_groups.tsv")
