library(tidyverse)
library(glue)
library(optparse)
# Renames headers and formats percolator output files for use with Ursgal's combine_pep_1_0_0.py script
# The ursgal script needs decoys as well as regular matches

read_psms <- function(file, is_decoy) {
  t <- read_tsv(file, col_names = FALSE) %>%
    mutate(X6 = unlist(lapply(X6, gsub, pattern = "\t", replacement = ","))) %>%
    `colnames<-`(.[1,]) %>%
    slice(-1) %>%
    select(c(score, `q-value`, posterior_error_prob, peptide)) %>%
    rename(PEP = posterior_error_prob) %>%
    mutate("Is_decoy" = is_decoy) %>%
    mutate(peptide = unlist(lapply(peptide, clean_peptide)))
  return(t)
}

read_tide <- function(file, is_decoy) {
  psm_header <- c("score", "q-value", "PEP", "peptide")
  t <- read_tsv(file) %>%
    select(c(
      "percolator score", `percolator q-value`,
      "percolator PEP", sequence
    )) %>%
    rename(all_of(c(
      score = "percolator score",
      `q-value` = "percolator q-value",
      PEP = "percolator PEP", peptide = "sequence"
    ))) %>%
    mutate("Is_decoy" = is_decoy) %>%
    mutate(peptide = unlist(lapply(peptide, clean_peptide)))
  return(t)
}

read_percolator <- function(valid_file, decoy_file, engine) {
  if (engine == "tide") {
    valid <- read_tide(valid_file, is_decoy = FALSE)
    decoy <- read_tide(decoy_file, is_decoy = TRUE)
    return(rbind(valid, decoy))
  } else {
    valid <- read_psms(valid_file, is_decoy = FALSE)
    decoy <- read_psms(decoy_file, is_decoy = TRUE)
    return(rbind(valid, decoy))
  }
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-m", "--matches"),
                       type = "character",
                       help = "File containing non-decoy matches only"
  )
  parser <- add_option(parser, c("-r", "--r_source"),
                       type = "character",
                       help = "Path to r scripts"
  )
  parser <- add_option(parser, c("-d", "--decoys"),
                       type = "character",
                       help = "File containing decoy matches only"
  )
  parser <- add_option(parser, c("-e", "--engine"),
                       help = "Set this if the files contain matches to
          proteins rather than peptides/psms"
  )
  parser <- add_option(parser, c("-o", "--output"),
                       type = "character",
                       help = "Output file name"
  )
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  combined <- read_percolator(args$matches, args$decoys, args$engine)
  write.csv(combined, file = args$output, row.names = FALSE)
}
