library(seqinr)
library(tidyverse)
library(glue)

# 1. Determine which proteins weren't matched by eggnog, then extract to a fasta file for further annotation by interpro
# 2. Merge identified eggnog proteins with previous metadata from percolator to identify their origin e.g. from which de novo search engine and obtain metadata

write_unmatched <- function(from_blast, eggnog_hits) {
  matched_by_eggnog <- (from_blast$ProteinId %in% eggnog_hits$ProteinId)
  not_matched_by_eggnog <- from_blast %>%
    filter(!(matched_by_eggnog))
  sequences <- as.list(not_matched_by_eggnog$seq)
  ids <- not_matched_by_eggnog$ProteinId
  ## write.fasta(sequences, ids, output_fasta)
  return(list(fasta_seqs = sequences, fasta_ids = ids, tsv = not_matched_by_eggnog))
}

main <- function(args) {
  unmatched_blast <- read_tsv(args$blast)
  anno_df <- read_tsv(args$annotations, skip = 4)
  seed_df <- read_tsv(args$seeds, skip = 5) %>% select(-c("evalue"))
  distinct_cols <- c(
    "#query", "evalue", "score", "bitscore",
    "pident", "qcov", "scov"
  )
  joined <- inner_join(anno_df, seed_df,
    by = join_by(x$`#query` == y$`#qseqid`)
  ) %>%
    select(-c("qstart", "qend", "sstart", "send", "sseqid")) %>%
    as_tibble()
  # Import the eggnog data and merge the two together, grouping
  # by queries that have all been matched to the same seed ortholog
  #   Don't do this (YET) because proteins/peptides matching to the same seed
  # ortholog doesn't necessarily mean that they belong to the same protein.
  # They may belong to separate proteins that are ALL homologs of the seed ortholog
  ## joined <- joined %>%
  ##   lapply(., as.character) %>%
  ##   as_tibble() %>%
  ##   group_by(seed_ortholog) %>%
  ##   mutate_at(., distinct_cols,
  ##     paste0,
  ##     collapse = ","
  ##   ) %>%
  ##   distinct() %>%
  ##   ungroup()
  with_blast <- inner_join(unmatched_blast, joined,
    by = join_by(x$ProteinId == y$`#query`)
  )
  ## transcriptome <- with_blast %>% filter(grepl("T"), ProteinId)
  transcriptome <- with_blast %>% filter(grepl("TRANSCRIPTOME"), header)
  # Get metadata for eggnog-identified proteins
  final <- mutate(with_blast,
    Anno_method = "eggNOG", GO = GOs,
    KEGG_ko = unlist(lapply(KEGG_ko, gsub,
      pattern = "ko:",
      replacement = ""
    ))
  ) %>% select(-GOs)
  anno_cols <- c(
    "seed_ortholog", "eggNOG_OGs",
    "max_annot_level", "COG_category", "Description",
    "Preferred_name", "GO", "EC", "KEGG_ko", "PFAMs",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
    "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction"
  )
  annotations <- final %>% select(any_of(c("ProteinId", "header", anno_cols)))
  meta <- final %>% select(-any_of(c(
    "header", anno_cols, "bitscore", "pident",
    "qcov", "scov", "score", "max_annot_lvl",
    "mass", "length"
  )))
  u <- write_unmatched(unmatched_blast, meta)
  return(list(meta_df = meta, anno_df = annotations, unmatched = u))
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-f", "--output_fasta"),
    type = "character",
    help = "Output fasta file name (unmatched eggnog proteins)"
  )
  parser <- add_option(parser, c("-u", "--output_unmatched"),
    type = "character",
    help = "Output unmatched tsv file name (unmatched eggnog proteins)"
  )
  parser <- add_option(parser, c("-m", "--output_meta"),
    type = "character",
    help = "Output metadata tsv file name"
  )
  parser <- add_option(parser, c("-o", "--output_anno"),
    type = "character",
    help = "Output annotation tsv file name"
  )
  parser <- add_option(parser, c("-s", "--seeds"),
    type = "character",
    help = "eggNOG seed file"
  )
  parser <- add_option(parser, c("-a", "--annotations"),
    type = "character",
    help = "eggNOG annotations file"
  )
  parser <- add_option(parser, c("-b", "--blast"),
    type = "character",
    help = "tsv containing unmatched blast hits originally given to eggnog"
  )
  ARGS <- parse_args(parser)
  m <- main(ARGS)
  write.fasta(m$unmatched$fasta_seqs, m$unmatched$fasta_ids, ARGS$output_fasta)
  write_delim(m$unmatched$tsv, ARGS$output_unmatched, delim = "\t")
  write_delim(m$meta_df, ARGS$output_meta, delim = "\t")
  write_delim(m$anno_df, ARGS$output_anno, delim = "\t")
}
