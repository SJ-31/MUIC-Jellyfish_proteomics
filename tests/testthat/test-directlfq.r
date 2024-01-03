library(glue)
source("./bin/directlfq_format2.r")
pth <- "./results/jellyfish/1-First_pass"
args <- list(
  path = glue("{pth}/Quantify/Mapped_scans"),
  mapping = "./results/jellyfish//msms_scans.tsv"
)

dlfq <- main(args)
