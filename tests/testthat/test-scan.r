library(testthat)
library(tidyverse)
library(glue)
bin <- "./bin"
source(glue("{bin}/get_scan_num.r"))

pth <- "./results/jellyfish/1-First_pass"
file_list <- list( # For tests
  comet = "../../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv",
  identipy = "../../results/test_manifest/1-First_pass/Percolator/identipy_percolator_psms.tsv",
  maxquant = "../../results/test_manifest/1-First_pass/MaxQuant/maxquant_all_pins.temp",
  msfragger = "../../results/test_manifest/1-First_pass/Percolator/msfragger_percolator_psms.tsv",
  metamorpheus = "../../results/test_manifest/1-First_pass/Metamorpheus/metamorpheus_AllPSMs.psmtsv",
  tide = "../../results/test_manifest/1-First_pass/Tide/tide-search.target.txt"
)

peps <- glue("{pth}/Quantify/Unmatched/unmatched_peptides.tsv") %>% read_tsv()
all_scans <- list.files(glue("{pth}/Quantify/Mapped_scans")) %>%
  lapply(., function(x) {
    paste0(glue("{pth}/Quantify/Mapped_scans/"), x)
  }) %>%
  unlist() %>%
  read_tsv()
