source("./bin/atleast2.r")
library(glue)
pth1 <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/"
pth <- "./tests/results/Percolator_proteins"
setwd(pth)
args <- list(
  seq_header_file = glue("{pth1}/Databases/seq-header_mappings.tsv"),
  pep_thresh = 1,
  fdr = 0.05
)
m <- main(args$seq_header_file, args$fdr, args$pep_thresh)
setwd("../../..")
