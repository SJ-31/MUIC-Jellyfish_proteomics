library(tidyverse)
library(ggplot2)
library(ggridges)
library(Peptides)
library(glue)

psm_file <- "./comet_percolator_psms.tsv"
prot_file <- "./comet_percolator_proteins.tsv"
map_file <- "../../../all_normal.fasta"
decoy_prefix <- "rev_"

psm_header <- c("PSMId", "score", "q-value", "posterior_error_prob", "peptide", "proteinIds")
engine_list <- c("comet", "msfragger", "metamorpheus", "identipy")
#TODO Add "maxquant" to the engine list

plot_mw_dist <- function(){
}

read_psms <- function(psm_file) {
  psms <- read.table(psm_file, sep = "\t", fill = TRUE, skip = 1,
                     header = FALSE) %>%
    select(c(1, 2, 3, 4, 5)) %>%
    as_tibble() %>%
    `colnames<-`(psm_header) %>%
    mutate(peptide = unlist(peptide %>% lapply(., gsub, pattern = "[^A-Z]*",
                                               replacement = ""))) %>%
    mutate(mw = unlist(lapply(peptide, mw))) %>%
    mutate(length = unlist(lapply(peptide, nchar)))
  return(psms)
}

get_tables <- function(path, is_psm) {
  if (is_psm == TRUE) {
    files <- paste0(path, paste0(engine_list, "_percolator_psms.tsv"))
    tables <- lapply(files, read_psms) %>% `names<-`(engine_list)
  } else {
    files <- paste0(path, paste0(engine_list, "_percolator_proteins.tsv"))
    tables <- lapply(files, read_psms) %>% `names<-`(engine_list)
  }
}

all_psms <- get_tables("./", TRUE)

all_psms <- seq_along(all_psms) %>% lapply(., function(x) {
  return(all_psms[[x]] %>% mutate(engine = all_psms[x]))
})

prot_header <- c("ProteinId", "ProteinGroupId", "q-value", "posterior_error_prob")

clear_decoys <- function(prot_list) {
  return(prot_list %>% gsub(glue(",*{decoy_prefix}.*($|,)"), "", .))
}
prot_table <- read.table(prot_file, sep = "\t", fill = TRUE, skip = 1,
                         header = FALSE) %>%
  select(c(1, 2, 3, 4)) %>%
  as_tibble() %>%
  `colnames<-`(prot_header) %>%
  mutate(ProteinId = unlist(ProteinId %>% lapply(., clear_decoys)))
