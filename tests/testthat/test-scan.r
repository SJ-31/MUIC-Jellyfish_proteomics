library(testthat)
library(tidyverse)
library(glue)
bin <- "./bin"
source(glue("{bin}/get_scan_num.r"))

pth <- "./results/jellyfish/1-First_pass"

peps <- glue("{pth}/Quantify/Unmatched/unmatched_peptides.tsv") %>% read_tsv()
all_scans <- list.files(glue("{pth}/Quantify/Mapped_scans")) %>%
  lapply(., function(x) {
    paste0(glue("{pth}/Quantify/Mapped_scans/"), x)
  }) %>%
  unlist() %>%
  read_tsv()

args <- list(
  msms_mapping = glue("./results/jellyfish/msms_scans.tsv"),
  input = glue("{pth}/Percolator/comet_percolator_psms.tsv"),
  engine = "comet",
  unmatched_peptides = glue("{pth}/Quantify/Unmatched/unmatched_peptides.tsv")
)
p <- read_engine_psms(args)
final <- merge_unmatched(p, args$unmatched_peptides)
