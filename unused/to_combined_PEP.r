#!/usr/bin/Rscript
library(tidyverse)
library(optparse)
## Renames headers and formats percolator output files for use with Ursgal's combine_pep_1_0_0.py script
## The ursgal script needs decoys as well as regular matches
parser <- OptionParser()
parser <- add_option(parser, c("-m", "--matches"), type="character",
                help="File containing non-decoy matches only")
parser <- add_option(parser, c("-d", "--decoys"), type="character",
                help="File containing decoy matches only")
parser <- add_option(parser, c("-p", "--protein_matches"), default = FALSE,
                     action = "store_true",
                     help = "Set this if the files contain matches to
        proteins rather than peptides/psms")
parser <- add_option(parser, c("-o", "--output"), type="character",
                     help = "Output file name")

clean_peptide <- function(modified_pep) {
  return(str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
         paste0(collapse = ""))
}

read_percolator <- function(filename, header, col_select, is_decoy, is_pep) {
  percolator_output <- read.table(filename, sep = "\t", fill = TRUE,
                                  skip = 1,
                                  header = FALSE) %>%
    as_tibble()
  percolator_output <- percolator_output %>%
    select(all_of(col_select)) %>%
    filter(!(grepl("[A-Za-z]", V2))) %>%
    filter(V2 != "") %>%
    `colnames<-`(header) %>%
    mutate("Is_decoy" = is_decoy)
  if (is_pep) {
  percolator_output <- mutate(percolator_output,
                              peptide = unlist(lapply(peptide, clean_peptide)))
   }
  return(percolator_output)
}

args <- parse_args(parser)
peptide_cols <- c("V2", "V3", "V4", "V5")
protein_cols <- c("V1", "V2", "V3", "V4")

# Testing purposes
## args <- list(matches = "../results/test_2combine/comet_percolator_psms.tsv",
##              decoys = "../results/test_2combine/comet_percolator_decoy_psms.tsv",
##              protein_matches = FALSE,
##              output = "../results/test_2combine/2combined.tsv")

psm_header <- c("score", "q-value", "PEP", "peptide")
percolator_protein_header <- c("ProteinId", "ProteinGroupID", "q-value", "PEP")

if (args$protein_matches) {
  valid <- read_percolator(args$matches, percolator_protein_header,
                           protein_cols, FALSE, FALSE)
  decoys <- read_percolator(args$decoys, percolator_protein_header,
                            protein_cols, TRUE, FALSE)
} else {
  valid <- read_percolator(args$matches, psm_header,
                           peptide_cols, FALSE, TRUE)
  decoys <- read_percolator(args$decoys, psm_header,
                            peptide_cols, TRUE, TRUE)
}

combined <- rbind(valid, decoys)
write.csv(combined, file = args$output, row.names = FALSE)
