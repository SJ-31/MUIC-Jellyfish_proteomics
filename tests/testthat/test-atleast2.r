library(testthat)
library(glue)
source("./bin/atleast2.r")
pth1 <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/"
pth <- "./tests/results/Percolator_proteins"
args <- list(
  seq_header_file = glue("{pth1}/Databases/seq-header_mappings.tsv"),
  path = pth,
  r_source = "./bin"
)
source(glue("{args$r_source}/helpers.r"))
output <- "~/atleast2_test.tsv"
m <- main(args$seq_header_file, args$path)
write_tsv(m, "./tests/testthat//output//intersected_searches.tsv")
