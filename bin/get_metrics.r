library(tidyverse)
library(ggplot2)
library(ggridges)
library(Peptides)
library(glue)


map_file <- "~/mapping.csv"
decoy_prefix <- "rev_"
file_path <- "../results/test_manifest/First_pass/Percolator/"

psm_header <- c("PSMId", "score", "q-value", "posterior_error_prob", "peptide", "proteinIds")
engine_list <- c("comet", "msfragger", "metamorpheus", "identipy")
prot_header <- c("ProteinId", "ProteinGroupId", "q.value", "posterior_error_prob")
# TODO Add "maxquant" to the engine list


read_psms <- function(psm_file) {
  psms <- read.table(psm_file,
    sep = "\t", fill = TRUE, skip = 1,
    header = FALSE
  ) %>%
    select(c(1, 2, 3, 4, 5)) %>%
    as_tibble() %>%
    `colnames<-`(psm_header)
  return(psms)
}

clear_decoys <- function(protein_list) {
  return(gsub(glue(",*{decoy_prefix}.*($|,)"), "", protein_list))
}

read_prot <- function(prot_file) {
  prot <- read.table(prot_file,
    sep = "\t", fill = TRUE, skip = 1,
    header = FALSE
  ) %>%
    as_tibble() %>%
    select(c(1, 2, 3, 4)) %>%
    `colnames<-`(prot_header) %>%
    mutate(ProteinId = unlist(ProteinId %>% lapply(., clear_decoys)))
  return(prot)
}

read_tide <- function(path, is_psm) {
  if (is_psm) {
    tide <- paste0(path, "tide_percolator_psms.tsv")
    tide_table <- read.table(tide,
      sep = "\t", fill = TRUE, skip = 1,
      header = FALSE
    ) %>%
      select(c(2, 7, 8, 9, 11)) %>%
      as_tibble() %>%
      `colnames<-`(psm_header)
    return(tide_table)
  }
}

format_psms <- function(psm_table) {
  psm_table <- psm_table %>%
    filter(!grepl(".*[A-Za-z]+.*", psm_table$score)) %>%
    mutate(peptide = unlist(peptide %>% lapply(., gsub,
      pattern = "[^A-Z]*",
      replacement = ""
    ))) %>%
    mutate(mw = unlist(lapply(peptide, mw))) %>%
    mutate(score = as.numeric(score)) %>%
    mutate(`q-value` = as.numeric(`q-value`)) %>%
    mutate(PSMId = as.character(PSMId)) %>%
    mutate(`posterior_error_prob` = as.numeric(`posterior_error_prob`)) %>%
    mutate(length = unlist(lapply(peptide, nchar)))
  return(psm_table)
}


split_duplicates <- function(dupe_table, index) {
  # Split a Percolator row containing duplicate protein ids into several rows, one for each id
  dupes <- dupe_table[index, ]$ProteinId %>% strsplit(",")
  split_dupes <- lapply(seq_along(dupes), function(x) {
    return(data.frame(
      ProteinId = dupes[x],
      ProteinGroupId = dupe_table[index, ]$ProteinGroupId,
      q.value = dupe_table[index, ]$q.value,
      posterior_error_prob = dupe_table[index, ]$posterior_error_prob
    ))
  })[[1]]
  colnames(split_dupes)[1] <- c("ProteinId")
  return(split_dupes)
}

sort_duplicates <- function(file_path) {
  # Read in a Percolator protein output file and sort duplicates
  table <- read_prot(file_path)
  duplicates <- table %>% filter(grepl(",", ProteinId))
  if (dim(duplicates)[1] == 0) {
    return(table)
  }
  table <- table %>% filter(!grepl(",", ProteinId))
  duplicates <- lapply(1:dim(duplicates)[1], split_duplicates, dupe_table = duplicates) %>%
    bind_rows()
  bound <- bind_rows(list(duplicates, table)) %>%
    as_tibble()
  return(bound)
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
    tables <- lapply(files, sort_duplicates) %>% `names<-`(engine_list)
    return(tables)
  }
}

protein_mappings <- read.delim(map_file, sep = "\t") %>% as_tibble()
prot_metrics <- function(protein_id, mode) {
  mapped <- protein_mappings %>% filter(id == protein_id)
  if (dim(mapped)[1] == 0) {
    return(NA)
  }
  if (mode == "length") {
    return(nchar(mapped$seq))
  } else if (mode == "mw") {
    return(mw(mapped$seq))
  }
}

## Completed
## all_psms <- get_tables(file_path, TRUE)
## tide_psms <- format_psms(read_tide(file_path, TRUE))
## all_psms$tide <- tide_psms
## all_psms <- seq_along(all_psms) %>% lapply(., function(x) {
##   return(all_psms[[x]] %>% mutate(engine = names(all_psms[x])))
## }) %>% bind_rows()

all_prot <- get_tables(file_path, FALSE)
all_prot <- seq_along(all_prot) %>%
  lapply(., function(x) {
    return(all_prot[[x]] %>% mutate(engine = names(all_prot[x])))
  }) %>%
  bind_rows() %>%
    mutate(mw = unlist(lapply(ProteinId, prot_metrics, mode = "mw"))) %>%
    mutate(length = unlist(lapply(ProteinId, prot_metrics, mode = "length")))

plot_mw_dist <- function(combined_table) {
  combined_table %>% ggplot(aes(x = mw, y = engine, fill = engine)) +
    ggridges::geom_density_ridges2()
}
