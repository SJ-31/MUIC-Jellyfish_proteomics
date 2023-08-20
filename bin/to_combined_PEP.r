library(tidyverse)
library(optparse)
## Renames headers and formats percolator output files for use with Ursgal's combine_pep_1_0_0.py script
## The ursgal script needs decoys as well as regular matches
args <- commandArgs(trailingOnly = TRUE)
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

read_percolator <- function(filename, header, col_select, is_decoy) {
  percolator_output <- read.table(filename, sep = "\t", fill = TRUE,
                                  skip = 1,
                                  header = FALSE) %>%
    as_tibble()
  percolator_output <- percolator_output %>%
    select(all_of(col_select)) %>%
    filter(!(grepl("[A-Za-z]", V2))) %>%
    filter(V2 != "") %>%
    `colnames<-`(header) %>%
    mutate("Is_decoy" = is_decoy) %>%
    mutate(peptide = unlist(lapply(peptide, gsub, pattern = "\\[.*\\]",
                                   replacement = ""))) %>%
    mutate(peptide = unlist(lapply(peptide, gsub, pattern = "\\.",
                                   replacement = ""))) %>%
    mutate(peptide = unlist(lapply(peptide, gsub, pattern = "-",
                                   replacement = "")))
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
                           protein_cols, FALSE)
  decoys <- read_percolator(args$decoys, percolator_protein_header,
                            protein_cols, TRUE)
} else {
  valid <- read_percolator(args$matches, psm_header,
                           peptide_cols, FALSE)
  decoys <- read_percolator(args$decoys, psm_header,
                            peptide_cols, TRUE)
}

combined <- rbind(valid, decoys)
write.csv(combined, file = args$output, row.names = FALSE)
