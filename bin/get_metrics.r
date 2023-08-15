library(tidyverse)
library(ggplot2)
library(ggridges)
library(Peptides)
library(glue)

map_file <- "../../../all_normal.fasta"
decoy_prefix <- "rev_"
file_path <- "./"

psm_header <- c("PSMId", "score", "q-value", "posterior_error_prob", "peptide", "proteinIds")
engine_list <- c("comet", "msfragger", "metamorpheus", "identipy")
#TODO Add "maxquant" to the engine list

plot_mw_dist <- function(combined_table){
  combined_table %>% ggplot(aes(x = mw, y = engine, fill = engine)) +
    ggridges::geom_density_ridges2()
}

read_psms <- function(psm_file) {
  psms <- read.table(psm_file, sep = "\t", fill = TRUE, skip = 1,
                     header = FALSE) %>%
    select(c(1, 2, 3, 4, 5)) %>%
    as_tibble() %>%
    `colnames<-`(psm_header)
  return(psms)
}

read_tide <- function(path, is_psm) {
  if (is_psm) {
    tide <- paste0(path, "tide_percolator_psms.tsv")
    tide_table <- read.table(tide, sep = "\t", fill = TRUE, skip = 1,
                             header = FALSE) %>%
      select(c(2, 7, 8, 9, 11)) %>%
      as_tibble() %>% `colnames<-`(psm_header)
    return(tide_table)
  }
}

format_psms <- function(psm_table) {
  psm_table <- psm_table %>%
    filter(!grepl(".*[A-Za-z]+.*", psm_table$score)) %>%
    mutate(peptide = unlist(peptide %>% lapply(., gsub, pattern = "[^A-Z]*",
                                               replacement = ""))) %>%
    mutate(mw = unlist(lapply(peptide, mw))) %>%
    mutate(score = as.numeric(score)) %>%
    mutate(`q-value` = as.numeric(`q-value`)) %>%
    mutate(PSMId = as.character(PSMId)) %>%
    mutate(`posterior_error_prob` = as.numeric(`posterior_error_prob`)) %>%
    mutate(length = unlist(lapply(peptide, nchar)))
  return(psm_table)
}

get_tables <- function(path, is_psm) {
  if (is_psm == TRUE) {
    files <- paste0(path, paste0(engine_list, "_percolator_psms.tsv"))
    tables <- lapply(files, read_psms) %>%
      lapply(., format_psms) %>%
      `names<-`(engine_list)
    return(tables)
  } else {
    files <- paste0(path, paste0(engine_list, "_percolator_proteins.tsv"))
    tables <- lapply(files, read_psms) %>% `names<-`(engine_list)
  }
}

all_psms <- get_tables(file_path, TRUE)
tide_psms <- format_psms(read_tide(file_path, TRUE))
all_psms$tide <- tide_psms
all_psms <- seq_along(all_psms) %>% lapply(., function(x) {
  return(all_psms[[x]] %>% mutate(engine = names(all_psms[x])))
}) %>% bind_rows()


prot_header <- c("ProteinId", "ProteinGroupId", "q-value", "posterior_error_prob")
# TODO: Finish up the graphing for the protein groups

clear_decoys <- function(prot_list) {
  return(prot_list %>% gsub(glue(",*{decoy_prefix}.*($|,)"), "", .))
}
prot_table <- read.table(prot_file, sep = "\t", fill = TRUE, skip = 1,
                         header = FALSE) %>%
  select(c(1, 2, 3, 4)) %>%
  as_tibble() %>%
  `colnames<-`(prot_header) %>%
  mutate(ProteinId = unlist(ProteinId %>% lapply(., clear_decoys)))
