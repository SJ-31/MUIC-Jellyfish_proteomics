library(testthat)
library(tidyverse)
library(glue)
bin <- "./bin"

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
  input = glue("{pth}/Percolator/msgf_percolator_psms.tsv"),
  engine = "metamorpheus",
  protein = glue("{pth}/Percolator/msgf_percolator_proteins.tsv"),
  unmatched_peptides = glue("{pth}/Quantify/Unmatched/unmatched_peptides.tsv")
)
source(glue("{bin}/get_scan_num.r"))

test_dlfq <- function() {
  args <- list(
    path = glue("{pth}/Quantify/Mapped_scans"),
    mapping = "./results/jellyfish/msms_scans.tsv"
  )
  source(glue("{bin}/directlfq_format2.r"))
  f <- main(args)
}

test_flfq <- function() {
  pth <- "./results/jellyfish/1-First_pass"
  args <- list(
    path = glue("{pth}/Quantify/Mapped_scans")
  )
  source(glue("{bin}/flashlfq_format2.r"))
  return(main(args))
}

f <- test_flfq()
