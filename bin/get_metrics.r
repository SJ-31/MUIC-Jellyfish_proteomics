library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)

read_tide <- function(tide_file) {
  tide <- read_tsv(tide_file) %>%
    select(c("percolator q-value", "percolator PEP", "sequence")) %>%
    dplyr::rename(PEP = `percolator PEP`) %>%
    dplyr::rename(q_value = `percolator q-value`) %>%
    dplyr::rename(peptide = sequence)
  return(tide)
}


read_psms <- function(file) {
  t <- read_tsv(file, col_names = FALSE) %>%
    mutate(X6 = unlist(lapply(X6, gsub, pattern = "\t", replacement = ","))) %>%
    `colnames<-`(.[1, ]) %>%
    slice(-1) %>%
    select(c(`q-value`, posterior_error_prob, peptide)) %>%
    dplyr::rename(PEP = posterior_error_prob) %>%
    dplyr::rename(q_value = `q-value`)
  return(t)
}

read_percolator_psms <- function(percolator_file) {
  engine <- gsub("_.*", "", percolator_file)
  percolator_file <- glue("{file_path}/{percolator_file}")
  if (engine == "tide") {
    p_tibble <- read_tide(percolator_file)
  } else {
    p_tibble <- read_psms(percolator_file)
  }
  p_tibble <- p_tibble %>%
    mutate(engine = engine) %>%
    mutate(peptide = unlist(lapply(peptide, cleanPeptide),
      use.names = FALSE
    )) %>%
    mutate(mw = mw(peptide)) %>%
    mutate(
      PEP = as.numeric(PEP),
      q_value = as.numeric(q_value)
    ) %>%
    as_tibble()
  return(p_tibble)
}

get_psms <- function(psm_list, pep_threshold, fdr) {
  all <- lapply(psm_list, read_percolator_psms) %>%
    bind_rows() %>%
    mutate(length = nchar(peptide)) %>%
    filter(PEP <= pep_threshold) %>%
    filter(q_value <= fdr)
}

clear_decoys <- function(protein_list) {
  decoy_prefix <- "rev_"
  return(gsub(glue(",*{decoy_prefix}.*($|,)"), "", protein_list))
}

read_prot <- function(prot_file) {
  engine <- gsub("_.*", "", prot_file)
  prot_file <- glue("{file_path}/{prot_file}")
  prot_header <- c("id", "q_value", "PEP")
  prot <- read_tsv(prot_file, col_names = FALSE) %>%
    as_tibble() %>%
    select(c(X1, X3, X4)) %>%
    slice(-1) %>%
    `colnames<-`(prot_header) %>%
    mutate(
      id = unlist(lapply(id, clear_decoys)),
      engine = engine,
      PEP = as.numeric(PEP),
      q_value = as.numeric(q_value)
    )
  return(prot)
}

split_duplicates <- function(dupe_table, index) {
  # Split a Percolator row containing duplicate protein ids into several rows, one for each id
  dupes <- dupe_table[index, ]$id %>%
    strsplit(",") %>%
    unlist(use.names = FALSE)
  others <- select(dupe_table[index, ], -id)
  return(tibble(id = dupes, others))
}

sort_duplicates <- function(file_path) {
  # Read in a Percolator protein output file and sort duplicates
  table <- read_prot(file_path)
  duplicates <- table %>% filter(grepl(",", id))
  if (dim(duplicates)[1] == 0) {
    return(table)
  }
  table <- table %>% filter(!grepl(",", id))
  duplicates <- lapply(1:dim(duplicates)[1], split_duplicates, dupe_table = duplicates) %>%
    bind_rows()
  bound <- bind_rows(list(duplicates, table)) %>%
    as_tibble() %>%
    mutate(PEP = as.numeric(PEP))
  return(bound)
}

get_prot <- function(prot_list, pep_threshold, fdr, mapping) {
  mapping <- read_tsv(mapping) %>%
    select(c(id, mass, length))
  all <- lapply(prot_list, sort_duplicates) %>%
    bind_rows() %>%
    filter(PEP <= pep_threshold) %>%
    filter(q_value <= fdr)
  joined <- inner_join(mapping, all, by = join_by(id)) %>%
    dplyr::rename(mw = mass)
  print(joined$mw)
  return(as_tibble(joined))
}

get_tables <- function(path, is_psm, pep_thresh, fdr, mapping) {
  if (is_psm == TRUE) {
    paths <- list.files(path, pattern = "*percolator_psms.tsv")
    return(get_psms(paths, pep_thresh, fdr))
  } else {
    paths <- list.files(path, pattern = "*percolator_proteins.tsv")
    return(get_prot(paths, pep_thresh, fdr, mapping))
  }
}
source(glue("{args$r_source}/helpers.r"))
