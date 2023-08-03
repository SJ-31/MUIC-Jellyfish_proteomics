library(tidyverse)
library(optparse)
## Renames headers and formats percolator output files for use with Ursgal's combine_pep_1_0_0.py script
args <- commandArgs(trailingOnly = TRUE)
## args <- c("comet_psms.tsv", "comet_decoy_psms.tsv", "combined_comet.csv")
## valid_file <- args[1]
## decoy_file <- args[2]
## combined_file <- args[3]
parser <- OptionParser()
parser <- add_option(parser, c("-m", "--matches"), type="character",
                help="File containing non-decoy matches only")
parser <- add_option(parser, c("-d", "--decoys"), type="character",
                help="File containing decoy matches only")
parser <- add_option(parser, c("-p", "--protein_matches"), default=FALSE,
                     action="store_true",
                     help="Set this if the files contain matches to
        proteins rather than peptides/psms")
parser <- add_option(parser, c("-q", "--maxquant"), type="character",
                help="Path to MaxQaunt msmsScans file, if applicable")
parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file name")
test_psms = c("--matches=comet_psms.tsv",
                                    "--decoys=comet_decoy_psms.tsv",
                                    "--output=combined_comet.csv"
                                    )
test_prot = c("--matches=comet_proteins.tsv",
                                    "--decoys=comet_decoy_proteins.tsv",
                                    "--output=combined_prot.csv", "-p"
                                    )
test_mq = c("--maxquant=combined_msms",
            "--output=maxquant_psms.csv"
                            )

read_percolator <- function(filename, header, col_select, is_decoy) {
  percolator_output <- read.table(filename, sep = "\t", fill = TRUE) %>%
    as_tibble() %>%
    select(all_of(col_select)) %>%
    filter(!(grepl("[A-Za-z]", V2))) %>%
    filter(V2 != "") %>%
    `colnames<-`(header)%>%
    mutate("Is decoy" = is_decoy)
  return(percolator_output)
}

## args <- parse_args(parser, args = test_mq)
args <- parse_args(parser)
peptide_cols <- c("V2", "V3", "V4", "V5")
protein_cols <- c("V1", "V2", "V3", "V4")
mq_cols <- c("V36", "V11", "V37", "V34", "V12")

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
