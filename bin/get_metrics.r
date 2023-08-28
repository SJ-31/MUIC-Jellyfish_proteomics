library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)

get_percolator_row <- function(row_index, percolator_lines,
                               names) {
  splits <- percolator_lines[row_index] %>%
    str_split("\t") %>%
    unlist()
  return(tibble(
    pep = splits[4],
    peptide = splits[5],
  ))
}

clean_peptide <- function(modified_pep) {
  clean_pep <- str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
    paste0(collapse = "")
  return(clean_pep)
}

read_tide <- function(tide_file) {
  tide <- read.delim(tide_file, sep = "\t") %>%
    select(c("percolator.PEP", "sequence")) %>%
    rename(pep = percolator.PEP) %>%
    rename(peptide = sequence)
  return(tide)
}

read_percolator_psms <- function(percolator_file) {
  engine <- gsub("_.*", "", percolator_file)
  percolator_file <- glue("{file_path}/{percolator_file}")
  if (engine == "tide") {
    p_tibble <- read_tide(percolator_file)
  } else {
  lines <- read_lines(percolator_file)
  p_tibble <- lapply(seq_along(lines)[-1], get_percolator_row,
    percolator_lines = lines
  ) %>%
    bind_rows()
  }
  p_tibble <- p_tibble %>%
    mutate(engine = engine) %>%
    mutate(peptide = unlist(lapply(peptide, clean_peptide), use.names = FALSE)) %>%
    mutate(mw = mw(peptide)) %>%
    mutate(pep = as.numeric(pep)) %>%
    as_tibble()
  return(p_tibble)
}

get_psms <- function(psm_list, pep_threshold) {
  all <- lapply(psm_list, read_percolator_psms) %>%
    bind_rows() %>%
    mutate(length = nchar(peptide)) %>%
    filter(pep <= pep_threshold)
}

clear_decoys <- function(protein_list) {
  decoy_prefix <- "rev_"
  return(gsub(glue(",*{decoy_prefix}.*($|,)"), "", protein_list))
}

read_prot <- function(prot_file) {
  engine <- gsub("_.*", "", prot_file)
  prot_file <- glue("{file_path}/{prot_file}")
  prot_header <- c("id", "pep")
  prot <- read.table(prot_file,
    sep = "\t", fill = TRUE, skip = 1,
    header = FALSE
  ) %>%
    as_tibble() %>%
    select(c(1, 4)) %>%
    `colnames<-`(prot_header) %>%
    mutate(id = unlist(id %>% lapply(., clear_decoys))) %>%
    mutate(engine = engine)
  return(prot)
}

split_duplicates <- function(dupe_table, index) {
  # Split a Percolator row containing duplicate protein ids into several rows, one for each id
  dupes <- dupe_table[index, ]$id %>% strsplit(",") %>% unlist(use.names = FALSE)
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
    mutate(pep = as.numeric(pep))
  return(bound)
}

get_prot <- function(prot_list, pep_threshold, mapping) {
  mapping <- read.delim(mapping, sep = "\t") %>%
    select(c(id, mass, length))
  all <- lapply(prot_list, sort_duplicates) %>%
    bind_rows() %>%
    filter(pep <= pep_threshold)
  joined <- inner_join(mapping, all, by = join_by(id)) %>%
                                       rename(mw = mass)
  return(as_tibble(joined))
}

get_tables <- function(path, is_psm, thresh, mapping) {
  if (is_psm == TRUE) {
    paths <- list.files(path, pattern = "*percolator_psms.tsv")
    return(get_psms(paths, thresh))
  } else {
    paths <- list.files(path, pattern = "*percolator_proteins.tsv")
    return(get_prot(paths, thresh, mapping))
  }
}
