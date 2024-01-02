library(tidyverse)
# Merge all percolator protein results into a single file
TARGET <- "ProteinId"
headers <- c(
  "ProteinId", "ProteinGroupId", "q.value", "posterior_error_prob",
  "peptideIds"
)

splitDuplicates <- function(dupe_table, index) {
  # Split a Percolator row containing duplicate protein ids into several rows, one for each id
  dupes <- dupe_table[index, ]$ProteinId %>%
    strsplit(",") %>%
    unlist(use.names = FALSE)
  others <- select(dupe_table[index, ], -ProteinId)
  return(tibble(ProteinId = dupes, others))
}

sort_duplicates <- function(file_path, fdr, pep_thresh) {
  # Read in a Percolator protein output file and sort duplicates, proteins
  # that are all known from the same set of peptides
  table <- read.delim(file_path, sep = "\t")
  engine <- gsub("_.*", "", file_path)
  duplicates <- table %>% filter(grepl(",", ProteinId))
  if (dim(duplicates)[1] == 0) {
    table <- table %>% mutate(ProteinGroupId = paste0(ProteinGroupId, engine))
    return(table)
  }
  table <- table %>% filter(!grepl(",", ProteinId))
  duplicates <- lapply(1:dim(duplicates)[1], splitDuplicates, dupe_table = duplicates) %>%
    bind_rows()
  bound <- bind_rows(list(duplicates, table)) %>%
    mutate(ProteinGroupId = paste0(ProteinGroupId, engine)) %>%
    as_tibble()
  bound <- bound %>% filter(`q.value` <= fdr)
  bound <- bound %>% filter(posterior_error_prob <= pep_thresh)
  return(bound)
}

get_matches <- function(file_name, target, fdr, pep_thresh) {
  engine <- gsub("_.*", "", file_name)
  results <- sort_duplicates(file_name, fdr, pep_thresh)
  matches <- results[[target]]
  engine_results <- list(matches)
  names(engine_results) <- engine
  return(engine_results)
}

intersect_engines <- function(files, map_file, fdr, pep_thresh) {
  mappings <- read.delim(map_file, sep = "\t")
  engi <- lapply(files, get_matches,
    target = TARGET,
    fdr = fdr,
    pep_thresh = pep_thresh
  ) %>%
    unlist(recursive = FALSE)
  tables <- lapply(files, sort_duplicates, fdr = fdr, pep_thresh = pep_thresh)
  combos <- combn(names(engi), 2)
  master_list <- lapply(1:ncol(combos), function(x) {
    # Obtain the intersection between every possible pair of engine searches
    engine_pair <- combos[, x]
    intersection <- intersect(
      engi[[engine_pair[1]]],
      engi[[engine_pair[2]]]
    )
    return(intersection)
  }) %>%
    unlist(use.names = FALSE) %>%
    unique()
  master_list <- master_list[!grepl("rev_", master_list)]
  matched_tables <- lapply(tables, function(x) {
    return(x[x[[TARGET]] %in% master_list, ])
  })
  mapped <- mappings[mappings$id %in% master_list, ]
  return(list("tables" = matched_tables, "map_list" = mapped))
}

merge_column <- function(column_name, dataframe) {
  # Combines information from all searches for a given match
  # E.g. for a given protein, the PEPs from all search engines will be formatted as a comma-delimited list
  # Will be useful to check the degree to which search engines agree with one
  # another
  # Will also help report which proteins a peptide has mapped to
  cols <- grep(column_name, colnames(dataframe))
  selected <- select(dataframe, all_of(cols))
  all_vals <- lapply(1:dim(selected)[1], function(x) {
    keep <- selected[x, ] %>% discard(~ all(is.na(.)))
    collapsed <- paste0(keep, collapse = ",")
    return(gsub(" ", "", collapsed))
  }) %>%
    unlist()
  return(all_vals)
}

main <- function(seq_header_file, output, fdr, pep_thresh) {
  files <- list.files(".", pattern = "*percolator.*") # This script will be run in a Nextflow process where
  # all the results files have been dumped
  # into the directory
  matched <- intersect_engines(files, seq_header_file, fdr, pep_thresh)
  matched_tables <- matched$tables
  mapped <- matched$map_list
  merged_tables <- reduce(matched_tables, full_join, by = TARGET)
  merged_tables <- lapply(headers, merge_column, dataframe = merged_tables) %>%
    `names<-`(headers) %>%
    as_tibble() %>%
    inner_join(mapped, by = join_by(x$ProteinId == y$id))
  write_delim(merged_tables, output, delim = "\t")
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-p", "--pep_thresh"),
    type = "character",
    help = "Posterior Error Probability threshold"
  )
  parser <- add_option(parser, c("-f", "--fdr"),
    type = "character",
    help = "FDR threshold"
  )
  parser <- add_option(parser, c("-m", "--seq_header_file"),
    type = "character",
    help = "Path to seq-header mapping"
  )
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  args <- parse_args(parser)
  main(args$seq_header_file, args$output, args$fdr, args$pep_thresh)
}
