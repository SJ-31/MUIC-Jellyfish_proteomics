library(tidyverse)
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-m", "--map_file"),
  type = "character",
  help = "Path to header mapping"
)
parser <- add_option(parser, c("-o", "--output"),
  type = "character",
  help = "Output file name"
)
args <- parse_args(parser)

files <- list.files(".", pattern = "*percolator.*") # This script will be run in a Nextflow process where
# all the results files have been dumped
# into the directory
map_file <- args$map_file
# Merge all percolator protein results into a single file
target <- "ProteinId"
headers <- c(
  "ProteinId", "ProteinGroupId", "q.value", "posterior_error_prob",
  "peptideIds"
)


split_duplicates <- function(dupe_table, index) {
  # Split a Percolator row containing duplicate protein ids into several rows, one for each id
  dupes <- dupe_table[index, ]$ProteinId %>% strsplit(",")
  split_dupes <- lapply(seq_along(dupes), function(x) {
    return(data.frame(
      ProteinId = dupes[x],
      ProteinGroupId = dupe_table[index, ]$ProteinGroupId,
      q.value = dupe_table[index, ]$q.value,
      posterior_error_prob = dupe_table[index, ]$posterior_error_prob,
      peptideIds = dupe_table[index, ]$peptideIds
    ))
  })[[1]]
  colnames(split_dupes)[1] <- c("ProteinId")
  return(split_dupes)
}
sort_duplicates <- function(file_path) {
  # Read in a Percolator protein output file and sort duplicates
  table <- read.delim(file_path, sep = "\t")
  engine <- gsub("_.*", "", file_path)
  duplicates <- table %>% filter(grepl(",", ProteinId))
  if (dim(duplicates)[1] == 0) {
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

filter_pep <- function(percolator_out, thresh) {
  filtered <- percolator_out[percolator_out$posterior_error_prob <= thresh, ]
  return(filtered)
}

get_matches <- function(file_name, target) {
  engine <- gsub("_.*", "", file_name)
  results <- sort_duplicates(file_name)
  ## results <- filter_pep(results, 0.1)
  matches <- results[[target]]
  engine_results <- list(matches)
  names(engine_results) <- engine
  return(engine_results)
}

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

matched_tables <- lapply(tables, function(x) {
  return(x[x[[target]] %in% master_list, ])
})

merge_column <- function(column_name, dataframe) {
  # Combines information from all searches for a given match
  # E.g. for a given protein, the PEPs from all search engines will be formatted as a comma-delimited list
  # Will be useful to check the degree to which search engines agree with one
  # another
  # Will also help report which proteins a peptide has mapped to
  cols <- grep(column_name, colnames(dataframe))
  selected <- select(dataframe, all_of(cols))
  all_vals <- lapply(1:dim(selected)[1], function(x) {
    collapsed <- paste0(selected[x, ], collapse = ",")
    return(gsub(" ", "", collapsed))
  }) %>%
    unlist()
  return(all_vals)
}

mappings <- read.delim(map_file, sep = "\t", col.names = c("id", "header"))
mapped <- mappings[mappings$id %in% master_list, ]
rm(mappings)

merged_tables <- reduce(matched_tables, full_join, by = target)
merged_tables <- lapply(headers, merge_column, dataframe = merged_tables) %>%
  `names<-`(headers) %>%
  as_tibble() %>%
  inner_join(mapped, by = join_by(x$ProteinId == y$id))

write_delim(merged_tables, args$output, delim = "\t")
mappings <- read.delim(map_file, sep = "\t", col.names = c("id", "header"))
mapped <- mappings[mappings$id %in% master_list, ]
rm(mappings)

merged_tables <- reduce(matched_tables, full_join, by = target)
merged_tables <- lapply(headers, merge_column, dataframe = merged_tables) %>%
  `names<-`(headers) %>%
  as_tibble() %>%
  inner_join(mapped, by = join_by(x$ProteinId == y$id))

write_delim(merged_tables, args$output, delim = "\t")
