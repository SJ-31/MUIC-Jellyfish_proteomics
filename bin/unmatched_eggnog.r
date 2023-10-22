library(seqinr)
library(tidyverse)
# Determine which proteins weren't matched by eggnog, then extract to a fasta file for further annotation by interpro

write_unmatched <- function(from_blast, eggnog_hits, output_fasta) {
  matched_by_eggnog <- (from_blast$ProteinId %in% eggnog_hits$`#query`)
  not_matched_by_eggnog <- from_blast %>%
    filter(!(matched_by_eggnog)) %>%
    distinct(seq, .keep_all = TRUE)
  sequences <- as.list(not_matched_by_eggnog$seq)
  ids <- not_matched_by_eggnog$ProteinId
  write.fasta(sequences, ids, output)
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-f", "--output_fasta"), type = "character",
                       help = "Output fasta file name")
  parser <- add_option(parser, c("-b", "--blast_input"), type = "character",
                       help = "TSV file of blast results given to eggnog")
  parser <- add_option(parser, c("-a", "--eggnog_annotations"), type = "character",
                       help = "emapper.annotations file")
  args <- parse_args(parser)
  write_unmatched(args$blast_input, args$eggnog_annotations, args$output_fasta)
}
