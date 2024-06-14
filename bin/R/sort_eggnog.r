library(seqinr)
library(tidyverse)
library(glue)


REPLACE_COMMA <- c("EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs", "GO")


# 1. Determine which proteins weren't matched by eggnog, then extract to a fasta file for further annotation by interpro
# 2. Merge identified eggnog proteins with previous metadata from percolator to identify their origin e.g. from which de novo search engine and obtain metadata
#

get_unmatched <- function(from_blast, eggnog_hits) {
  matched_by_eggnog <- (from_blast$ProteinId %in% eggnog_hits$ProteinId)
  not_matched_by_eggnog <- from_blast %>%
    filter(!(matched_by_eggnog))
  sequences <- as.list(not_matched_by_eggnog$seq)
  ids <- not_matched_by_eggnog$ProteinId
  if (nrow(not_matched_by_eggnog) == 0) {
    return(NULL)
  }
  return(list(fasta_seqs = sequences, fasta_ids = ids, tsv = not_matched_by_eggnog))
}

#' Find the true start of an eggnog annotations
#' file
find_start <- function(file) {
  lines <- readLines(file, n = 50)
  counter <- 0
  for (i in seq_len(length(lines))) {
    split <- str_split_1(lines[i], "\t")
    if (split[1] == "#qseqid" || split[1] == "#query") {
      break
    }
    counter <- counter + 1
  }
  return(counter)
}


main <- function(args) {
  unmatched_blast <- read_tsv(args$blast)
  anno_df <- read_tsv(args$annotations, skip = find_start(args$annotations))
  seed_df <- read_tsv(args$seeds, skip = find_start(args$seeds)) %>% select(-evalue)
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
  with_blast <- inner_join(unmatched_blast, joined,
    by = join_by(x$ProteinId == y$`#query`)
  )
  # Get metadata for eggnog-identified proteins
  final <- mutate(with_blast,
    inferred_by = "eggNOG",
    GO = GOs,
    KEGG_ko = unlist(lapply(KEGG_ko, gsub,
      pattern = "ko:",
      replacement = "",
    ))
  ) %>% select(-c(
    "GOs", "bitscore", "pident", "qcov", "scov",
    "score", "evalue", "bitscore", "max_annot_lvl"
  ))
  for (to_replace in REPLACE_COMMA) {
    replaced <- purrr::map_chr(
      final[[to_replace]],
      \(x) str_replace_all(x, ",", ";")
    )
    final <- dplyr::mutate(final, !!to_replace := replaced)
  }
  u <- get_unmatched(unmatched_blast, final)
  return(list(all = final, unmatched = u))
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
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "output file name"
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
  if (!is.null(m$unmatched)) {
    write.fasta(m$unmatched$fasta_seqs, m$unmatched$fasta_ids, ARGS$output_fasta)
    write_delim(m$unmatched$tsv, ARGS$output_unmatched, delim = "\t")
  } else {
    cat(file = "annotation_complete.fasta")
  }
  write_delim(m$all, ARGS$output, delim = "\t")
}
