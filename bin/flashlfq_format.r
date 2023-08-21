library(tidyverse)
library(Peptides)
library(optparse)
library(glue)

## parser <- OptionParser()
## parser <- add_option(parser, c("--metamorpheus"), type="character",
##                 help="metamorpheus psms")
## parser <- add_option(parser, c("--identipy"), type="character",
##                 help="identipy psms")
## parser <- add_option(parser, c("--maxquant"), type="character",
##                      help = "maxquant psms")
## parser <- add_option(parser, c("--msfragger"), type="character",
##                      help = "msfragger psms")
## parser <- add_option(parser, c("--comet"), type="character",
##                      help = "comet psms")
## parser <- add_option(parser, c("--tide"), type="character",
##                      help = "tide psms")
## parser <- add_option(parser, c("-o", "--output"), type="character",
##                      help = "Output file name")
## parser <- add_option(parser, c("-m", "--msms_mapping"), type="character",
##                      help = "MsMs mapping file")
## args <- parse_args(parser)
##

comet_input <- "../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv"
msms_ref <- "../results/CiCs_metrics.tsv"
# One of msfragger, comet, maxquant, identipy, metamorpheus, tide

group_prot <- function(prot_vec) {
  grouped <- paste0(prot_vec, collapse = ";")
  split_groups <- grouped %>%
    strsplit(";") %>%
    unlist(use.names = FALSE) %>%
    unique()
  return(paste0(split_groups, collapse = ";"))
}

get_percolator_row <- function(row_index, percolator_lines) {
    splits <- percolator_lines[row_index] %>% str_split("\t") %>% unlist()
    return(tibble(PSMId = splits[1],
            peptide = splits[5],
            `Protein Accession` = group_prot(splits[6: length(splits)])))
}

read_percolator <- function(percolator_file) {
  lines <- read_lines(percolator_file)
  p_tibble <- lapply(seq_along(lines)[-1], get_percolator_row,
                       percolator_lines = lines) %>% bind_rows()
  return(p_tibble)
}

##  Functions for formatting scan number
##    The scan number will let you match retension time
comet_scans <- function(comet_id_str) {
  path_removed <- gsub(".*/", "", comet_id_str)
  name_cleared <- strsplit(path_removed, "_") %>% unlist(use.names = FALSE)
  return(glue("{name_cleared[1]}.{name_cleared[2]}"))
}

maxquant_scans <- function(mq_id_str) {
  match_result <- regexec("(.*\\..*)\\..*", mq_id_str)
  groups <- regmatches(mq_id_str, match_result)[[1]]
  return(groups[2])
}

msfragger_scans <- function(msfragger_id_str) {
  match_result <- regexec("(.*\\..*)\\..*\\..*", msfragger_id_str)
  groups <- regmatches(msfragger_id_str, match_result)[[1]]
  return(groups[2])
}


get_file_name <- function(scan) {
  return(gsub("\\..*", "", scan))
}

clean_peptide <- function(modified_pep) {
  return(str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
         paste0(collapse = ""))
}

mapping <- read.delim(msms_ref, sep = "\t")

read_engine_psms <- function(percolator_input, engine) {
  # Read the percolator psm file from <engine> and format the
  #     the output according to flashlfq
  #     This relies on all functions defined above
  psms <- read_percolator(percolator_input)
  if ( engine == "comet" ) {
  psms <- psms %>% mutate(scan = unlist(lapply(PSMId, comet_scans)))
  } else if ( engine == "maxquant" ) {
  psms <- psms %>% mutate(scan = unlist(lapply(PSMId, maxquant_scans)))
  } else if ( engine == "msfragger" ) {
  psms <- psms %>% mutate(scan = unlist(lapply(PSMId, msfragger_scans)))
  } else if ( engine == "tide" ) {

  } else if ( engine == "metamorpheus" ) {

  }
  ## Nothing to do for identipy
  ## TODO: Add formmatting for metamorpheus
  ## TODO: Add formmatting for tide
  psms <- psms %>%
  mutate(`File Name` = unlist(lapply(scan, get_file_name))) %>%
  mutate(`Full Sequence` = peptide) %>%
  mutate(`Base Sequence` = unlist(lapply(peptide, clean_peptide)))%>%
  mutate(`Peptide Monoisotopic Mass` = unlist(
          lapply(`Base Sequence`, mw, monoisotopic = TRUE))) %>%
    select(c("File Name", "scan", "Base Sequence", "Full Sequence",
             "Peptide Monoisotopic Mass",
             "Protein Accession"))
  new_cols <- c("Scan Retension Time", "Precursor Charge")
  old_cols <- c("retensionTime", "precursorCharge")
  psms <- inner_join(psms, mapping, by = join_by(scan == scanNum)) %>%
    rename_with( ~new_cols, all_of(old_cols)) %>%
    select(c("File Name", "Scan Retension Time", "Precursor Charge",
             "Base Sequence", "Full Sequence", "Peptide Monoisotopic Mass",
             "Protein Accession"))
  return(psms)
}

## file_list <- list(comet = args$comet, identipy = args$identipy,
##                   tide = args$tide, metamorpheus = args$metamorpheus,
##                   msfragger = args$msfragger, maxquant = args$maxquant)

file_list <- list(comet = "../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv",
identipy = "../results/test_manifest/1-First_pass/Percolator/identipy_percolator_psms.tsv", maxquant =  "../results/test_manifest/1-First_pass/Percolator/maxquant_percolator_psms.tsv", msfragger = "../results/test_manifest/1-First_pass/Percolator/msfragger_percolator_psms.tsv")

unique_peps <- vector()
comet <- read_engine_psms(comet_input, "comet")
unique_peps <- append(unique_peps, comet$`Base Sequence`)
