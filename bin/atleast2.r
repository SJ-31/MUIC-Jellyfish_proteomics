library(tidyverse)
library(optparse)
# Merge all percolator protein results into a single file
target <- "ProteinId"
headers <- c(
  "ProteinId", "ProteinGroupId", "q.value", "posterior_error_prob",
  "peptideIds"
)

split_duplicates <- function(dupe_table, index) {
  # Split a Percolator row containing duplicate protein ids into several rows, one for each id
  dupes <- dupe_table[index, ]$ProteinId %>%
    strsplit(",") %>%
    unlist(use.names = FALSE)
  others <- select(dupe_table[index, ], -ProteinId)
  return(tibble(ProteinId = dupes, others))
}

sort_duplicates <- function(file_path) {
  # Read in a Percolator protein output file and sort duplicates
  table <- read.delim(file_path, sep = "\t")
  engine <- gsub("_.*", "", file_path)
  duplicates <- table %>% filter(grepl(",", ProteinId))
  if (dim(duplicates)[1] == 0) {
    table <- table %>% mutate(ProteinGroupId = paste0(ProteinGroupId, engine))
    return(table)
  }
  table <- table %>% filter(!grepl(",", ProteinId))
  duplicates <- lapply(1:dim(duplicates)[1], split_duplicates, dupe_table = duplicates) %>%
    bind_rows()
  bound <- bind_rows(list(duplicates, table)) %>%
    mutate(ProteinGroupId = paste0(ProteinGroupId, engine)) %>%
    as_tibble()
  return(bound)
}

get_matches <- function(file_name, target) {
  engine <- gsub("_.*", "", file_name)
  results <- sort_duplicates(file_name)
  matches <- results[[target]]
  engine_results <- list(matches)
  names(engine_results) <- engine
  return(engine_results)
}

intersect_engines <- function(files, map_file) {
  mappings <- read.delim(map_file, sep = "\t")
  engi <- lapply(files, get_matches, target = target) %>%
    unlist(recursive = FALSE)
  tables <- lapply(files, sort_duplicates)
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
    return(x[x[[target]] %in% master_list, ])
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
    keep <- selected[x, ] %>% discard(~all(is.na(.)))
    collapsed <- paste0(keep, collapse = ",")
    return(gsub(" ", "", collapsed))
  }) %>%
    unlist()
  return(all_vals)
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-m", "--seq_header_file"),
    type = "character",
    help = "Path to seq-header mapping"
  )
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  args <- parse_args(parser)
  ## args <- list(map_file = "../workflow/results/ND_jellyfish/Databases/header_mappings.tsv", output = "~/intersect_test.tsv")
  files <- list.files(".", pattern = "*percolator.*") # This script will be run in a Nextflow process where
  # all the results files have been dumped
  # into the directory
  matched <- intersect_engines(files, args$seq_header_file)
  matched_tables <- matched$tables
  mapped <- matched$map_list
  merged_tables <- reduce(matched_tables, full_join, by = target)
  merged_tables <- lapply(headers, merge_column, dataframe = merged_tables) %>%
    `names<-`(headers) %>%
    as_tibble() %>%
    inner_join(mapped, by = join_by(x$ProteinId == y$id))
  write_delim(merged_tables, args$output, delim = "\t")
}
