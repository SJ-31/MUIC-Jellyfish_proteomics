library(tidyverse)
library(glue)
pth <- "./results/jellyfish/1-First_pass/Quantify/Mapped_scans/"

Mode <- function(x) {
  counts <- table(x)
  max <- max(counts)
  possibilities <- keep(counts, \(x) x == max) %>% names()
  if (length(possibilities) > 1) {
    return(sample(possibilities, 1))
  }
  return(possibilities)
}

## resolveScan <- function(grouped_df) {
##   p <- pick(base_peptide)$base_peptide
##   up <- unique(p)
##   tab <- tabulate(match(p, up))
##   abund <- up[tab == max(tab)]
##   current <- grouped_df[cur_group_rows(), ]
##   choices <- filter(current, base_peptide %in% abund)
##   if (nrow(choices) > 1) {
##     return(sample_n(choices, 1))
##   }
##   return(choices)
## }

mapped_scans <- list.files(pth) %>%
  lapply(., function(x) {
    read_tsv(glue("{pth}/{x}"))
  }) %>%
  bind_rows() %>%
  filter(!is.na(protein))

## proteinids <- read_lines("./tests/testthat/output//all_proteinids.txt")
