library(tidyverse)
library(testthat)
options(Debug = TRUE)
source("./bin/sort_eggnog.r")
pth <- "./results/jellyfish//1-First_pass/Unmatched"
args <- list(
  blast = glue("{pth}/BLAST/jellyfish_blast_unmatched.tsv"),
  peptides = glue("{pth}/unmatched_peptides.tsv"),
  annotations = glue("{pth}/eggNOG/jellyfish.emapper.annotations"),
  seeds = glue("{pth}/eggNOG/jellyfish.emapper.seed_orthologs")
)
m <- main(args)
