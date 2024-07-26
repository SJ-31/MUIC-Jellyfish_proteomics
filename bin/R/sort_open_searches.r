library(tidyverse)
library(glue)
#'
#' Sort open searches, including combining different peptides
#' of the same protein together
#'

dropNA <- function(vector) {
  return(vector[!is.na(vector)])
}

cleanNA <- function(vector) {
  cleaned <- lapply(vector, function(x) {
    if (is.na(x)) {
      return("")
    } else {
      return(x)
    }
  })
  return(unlist(cleaned, use.names = FALSE))
}

cleanUp <- function(path) {
  # Format protein rows with multiple peptides (duplicates)
  # Separate duplicates into separate rows
  engine <- gsub(".*/", "", path) %>% gsub("_.*", "", .)
  df <- read_tsv(path) %>%
    mutate(ProteinGroupId = paste0(ProteinGroupId, engine))

  df <- df |> mutate(peptideIds = fill_peptide_gaps(peptideIds))
  df <- separate_longer_delim(df, "ProteinId", ",") # Percolator uses comma to delimit protein ids
  return(df)
}


main <- function(args) {
  mapped <- read_tsv(args$seq_header_file)
  files <- list.files(args$path,
    pattern = "*percolator_proteins.tsv",
    full.names = TRUE
  )
  cleaned <- files %>% lapply(., cleanUp)
  grouped <- bind_rows(cleaned) %>%
    group_by(ProteinId) %>%
    mutate_at(., vars(-group_cols()), paste0, collapse = ";") %>%
    distinct() %>%
    ungroup() %>%
    inner_join(mapped, by = join_by(x$ProteinId == y$id))
  return(grouped)
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library(tidyverse)
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-r", "--r_source"),
    type = "character"
  )
  parser <- add_option(parser, c("-p", "--path"),
    type = "character",
    help = "Path to Protein tsvs"
  )
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  parser <- add_option(parser, c("-m", "--seq_header_file"),
    type = "character",
    help = "Path to seq-header mapping"
  )
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  m <- main(args) %>% mutate(ID_method = "open")
  write_delim(m, args$output, delim = "\t")
}
