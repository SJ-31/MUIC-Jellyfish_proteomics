#'
#' Sort open searches, including combining different peptides
# of the same protein together
#'

joinMods <- function(peptides) {
  pattern <- "([A-Z\\]]) ([A-Z\\[])"
  splits <- str_replace_all(peptides, pattern, "\\1,\\2") %>%
    unlist(use.names = FALSE)
  return(paste0(splits, collapse = ","))
}

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

sortDuplicates <- function(df) {
  dupes <- df %>% filter(grepl(",", df$ProteinId))
  if (dim(dupes)[1] == 0) {
    return(df)
  }
  no_dupes <- df[!(df$ProteinId %in% dupes$ProteinId), ]
  split <- lapply(seq_len(dim(dupes)[1]), splitDuplicates,
    dupe_table = dupes
  ) %>% bind_rows()
  return(bind_rows(split, no_dupes))
}

cleanUp <- function(path) {
  # Format protein rows with multiple peptides (duplicates)
  # Separate duplicates into separate rows
  # Filter by q-value
  engine <- gsub("_.*", "", path)
  df <- read_tsv(path) %>%
    mutate(ProteinGroupId = paste0(ProteinGroupId, engine))
  df <- mutate(df, peptideIds = unlist(lapply(df$peptideIds, joinMods),
    use.names = FALSE
  ))
  print(colnames(df))
  df <- sortDuplicates(df)
  return(df)
}


main <- function(args) {
  oldwd <- getwd()
  setwd(args$r_source)
  source("./atleast2.r")
  setwd(oldwd)
  mapped <- read_tsv(args$seq_header_file)
  files <- list.files(oldwd,
    pattern = "*percolator_proteins.tsv"
  )
  cleaned <- files %>% lapply(., cleanUp)
  grouped <- bind_rows(cleaned) %>%
    group_by(ProteinId) %>%
    mutate_at(., vars(-group_cols()), paste0, collapse = ",") %>%
    distinct() %>%
    ungroup() %>%
    inner_join(mapped, by = join_by(x$ProteinId == y$id))
  write_delim(grouped, args$output, delim = "\t")
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library(tidyverse)
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-r", "--r_source"),
    type = "character"
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
  main(args)
}
