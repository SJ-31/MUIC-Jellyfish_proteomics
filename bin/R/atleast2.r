library(glue)
library(tidyverse)

# Merge all percolator protein results into a single file
TARGET <- "ProteinId"
headers <- c(
  "ProteinId", "ProteinGroupId", "q.value", "posterior_error_prob",
  "peptideIds"
)

ID_LIST <- list()

get_matches <- function(file_name) {
  engine <- gsub(".*/", "", file_name) %>% gsub("_.*", "", .)
  results <- tibbleDuplicateAt(read_tsv(file_name), "ProteinId", ",") %>%
    mutate(ProteinGroupId = paste0(ProteinGroupId, engine))
  matches <- results[[TARGET]]
  engine_results <- list(matches)
  names(engine_results) <- engine
  ID_LIST <<- c(ID_LIST, engine_results)
  return(results)
}

intersect_engines <- function(files, map_file) {
  mappings <- read.delim(map_file, sep = "\t")
  tables <- lapply(files, get_matches)
  combos <- combn(names(ID_LIST), 2)
  master_list <- lapply(1:ncol(combos), function(x) {
    # Obtain the intersection between every possible pair of engine searches
    engine_pair <- combos[, x]
    intersection <- intersect(
      ID_LIST[[engine_pair[1]]],
      ID_LIST[[engine_pair[2]]]
    )
    return(intersection)
  }) %>%
    unlist(use.names = FALSE) %>%
    unique()
  master_list <- master_list[!grepl("rev_", master_list)]
  matched_tables <- lapply(tables, function(x) {
    return(x[x[[TARGET]] %in% master_list,])
  })
  mapped <- mappings[mappings$id %in% master_list,]
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
    keep <- selected[x,] %>% discard(~all(is.na(.)))
    collapsed <- paste0(keep, collapse = ";")
    return(gsub(" ", "", collapsed))
  }) %>%
    unlist()
  return(all_vals)
}

main <- function(seq_header_file, path) {
  # This script will be run in a Nextflow process where
  # all the results files have been dumped
  # into the directory
  files <- list.files(path, pattern = "*percolator.*", full.names = TRUE)
  matched <- intersect_engines(files, seq_header_file)
  matched_tables <- matched$tables
  mapped <- matched$map_list
  merged_tables <- reduce(matched_tables, full_join, by = TARGET)
  merged_tables <- lapply(headers, merge_column, dataframe = merged_tables) %>%
    `names<-`(headers) %>%
    as_tibble() %>%
    inner_join(mapped, by = join_by(x$ProteinId == y$id))
  return(merged_tables)
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-m", "--seq_header_file"),
                       type = "character",
                       help = "Path to seq-header mapping"
  )
  parser <- add_option(parser, c("-p", "--path"),
                       type = "character",
                       help = "Path to Protein tsvs"
  )
  parser <- add_option(parser, c("-r", "--r_source"),
                       type = "character",
                       help = "R source directory"
  )
  parser <- add_option(parser, c("-o", "--output"),
                       type = "character",
                       help = "Output file name"
  )
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  m <- main(args$seq_header_file, args$path)
  write_delim(m, args$output, delim = "\t")
}
