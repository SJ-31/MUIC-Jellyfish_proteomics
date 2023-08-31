library(MSnbase)
library(tidyverse)
library(optparse)
library(glue)

get_percolator_row <- function(row_index, percolator_lines) {
  splits <- percolator_lines[row_index] %>%
    str_split("\t") %>%
    unlist()
  return(tibble(
    pep = splits[4],
    peptide = splits[5],
  ))
}


read_percolator <- function(percolator_file) {
  lines <- read_lines(percolator_file)
  p_tibble <- lapply(seq_along(lines)[-1], get_percolator_row,
    percolator_lines = lines
  ) %>%
    bind_rows()
  return(p_tibble)
}

clean_peptide <- function(modified_pep) {
  clean_pep <- str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
    paste0(collapse = "")
  return(clean_pep)
}

join_scans_with_pep <- function(scan_file, perc_file, engine) {
  mapped_scans <- read.delim(scan_file, sep = "\t") %>% as_tibble()
  if (engine == "tide") {
    perc <- read.delim(perc_file, sep = "\t") %>%
      select(c(percolator.PEP, sequence)) %>%
      rename(pep = percolator.PEP) %>%
      rename(peptide = sequence)
  } else {
    perc <- read_percolator(perc_file) %>%
      mutate(peptide = unlist(lapply(peptide, clean_peptide)))
  }
  joined <- inner_join(mapped_scans, perc,
    by = join_by(base_peptide == peptide)
  ) %>%
    select(c(file, scan, pep)) %>%
    mutate(scan = unlist(gsub(".*\\.", "", scan)))
  return(joined)
}

spectra_names <- function(df) {
  # Format scan names to match the default scan names of the Spectra class
  df <- df %>%
    mutate(scan = unlist(lapply(scan, function(x) {
      return(glue("F1.S{x}"))
    })))
  return(df)
}


filter_scans <- function(psm_scans, mapping, pep_thresh) {
  psm_scans <- psm_scans %>% filter(pep <= pep_thresh)
  matched_psms <- psm_scans$scan
  unmatched <- psm_scans %>%
    filter(!(scan %in% mapping)) %>%
    select(c(file, scan))
  unmatched <- spectra_names(unmatched)
  return(unmatched)
}


filter_msms <- function(file_name, msms_path, df) {
  df <- filter(df, file == file_name)
  mzml <- readMSData(msms_path)
  sp <- spectra(mzml)
  unmatched_mzml <- mzml[names(sp) %in% df$scan]
  rm(mzml)
  return(unmatched_mzml)
}

filter_from_mzML <- function(msms_files, scan_df) {
  scans_left <- lapply(names(msms_files), function(x) {
    filter_msms(x, msms_files[[x]], scan_df)
  }) %>% `names<-`(names(msms_files))
  lapply(names(scans_left), function(x) {
    writeMSData(scans_left[[x]], file = glue("filtered_{x}.mzML"))
  })
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  parser <- OptionParser()
  parser <- add_option(parser, c("-m", "--ms_mapping"))
  parser <- add_option(parser, c("-s", "--scan_file"))
  parser <- add_option(parser, c("-r", "--percolator_file"))
  parser <- add_option(parser, c("-e", "--engine"))
  parser <- add_option(parser, c("-p", "--pep_thresh"))
  parser <- add_option(parser, c("-z", "--mzML_path"))
  args <- parse_args(parser)
  ms_map <- read.delim(args$ms_mapping, sep = "\t") %>% as_tibble()
  joined <- join_scans_with_pep(args$scan_file, args$percolator_file, args$engine)
  filtered <- filter_scans(joined, ms_map, as.numeric(args$pep_thresh))
  path <- args$mzML_path
  msms_files <- as.list(list.files(path,
    pattern = ".*\\.mzML"
  )) %>%
    `names<-`(gsub("\\..*", "", unlist(.)))
  filter_from_mzML(msms_files, filtered)
}
