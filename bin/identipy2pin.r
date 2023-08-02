library(tidyverse)
files <- commandArgs(trailingOnly = TRUE)
input <- files[1]
pin_file <- files[2]
sample_id <- files[3]

identipy <- read.delim(input, sep = "\t") %>% as_tibble()
new_names <- c(Peptide = "Sequence")
find_scan <- function(identipy_title) {
  match <- regexec(pattern = "scan=(.*)", identipy_title)
  scan <- regmatches(identipy_title, match)
  return(as.integer(scan[[1]][2]))
}
decoy_label <- function(identipy_proteins) {
  split <- unlist(strsplit(identipy_proteins, ";", fixed = TRUE))
  if (substr(split[1], 1, 4) == "rev_") {
    return(-1)
  }
  return(1)
}

split_to_tab <- function(identipy_proteins) {
  to_tab <- unlist(strsplit(identipy_proteins, ";", fixed = TRUE))
  return(paste0(to_tab, collapse = "\t"))
}

identipy <- identipy %>%
  mutate(Proteins = unlist(lapply(Proteins, gsub,
    pattern = ";",
    replacement = "\t"
  ))) %>%
  mutate(ScanNr = unlist(lapply(Title, find_scan))) %>%
  mutate(Label = unlist(lapply(Proteins, decoy_label))) %>%
  mutate(ScanID = unlist(lapply(ScanNr, function(x) {
    return(paste0(sample_id, ".", x))
  }))) %>%
  rename(., all_of(new_names)) %>%
  select(-c(Title, Modified.sequence))

identipy <- identipy[, c(
  19, 18, 17, 16, 15, 14, 13, 11, 9,
  8, 7, 6, 5, 4, 3, 2, 1, 12, 10
)]



pin <- seq_along(identipy$ScanNr) %>%
  lapply(., function(x) {
    return(identipy[x, ] %>% paste0(collapse = "\t"))
  }) %>%
  unlist()

cat(c(colnames(identipy), "\n"), file = pin_file, sep = "\t")
write_lines(pin, pin_file, append = TRUE)
