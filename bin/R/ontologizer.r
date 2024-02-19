library(tidyverse)

goVector <- function(df, column, filter) {
  # Return a flattened vector of GO terms from df
  df <- df %>% filter(!!as.symbol(column) == filter)
  gos <- df$GO %>%
    lapply(., strsplit, split = "\\||;|,") %>%
    unlist() %>%
    discard(is.na) %>%
    discard(!grepl("GO:", .))
}

writeOz <- function(file, df) {
  # Write ontologizer-formatted df to file
  header <- "GoStat IDs Format Version 1.0"
  write_lines(df, file)
  write_tsv(df, file, append = TRUE)
}

ozFormat <- function(df, with_filter, col, val) {
  # Format df as input for ontologizer program
  # Optional filtering by specific columns and values
  if (with_filter) {
    df <- df %>% filter(!!as.symbol(col) == val)
  }
  prepped <- df %>%
    select(c(ProteinId, GO)) %>%
    filter(!is.na(GO)) %>%
    mutate(GO = sapply(GO, function(x) {
      s <- unlist(strsplit(x, "\\||;|,")) %>% discard(!grepl("GO:", .))
      x <- sapply(s, gsub,
        pattern = "_.*", replacement = "",
        USE.NAMES = FALSE
      )
      return(paste0(x, collapse = ",")) # comma delimiter is required for ontologizer formatting
    }, USE.NAMES = FALSE)) %>%
    filter(GO != "NA")
  return(prepped)
}

prep <- function(args) {
  combined <- read_tsv(args$input)
  u <- ozFormat(combined, FALSE)
  io <- ozFormat(combined, TRUE, "ID_method", "open")
  sa <- ozFormat(filter(
    combined, Anno_method == "blast" | Anno_method == "intial_database",
  ), FALSE)
  return(list(
    universe = u,
    id_open = io,
    standard_annotation = sa
  ))
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"))
  args <- parse_args(parser)
  m <- prep(args)
  writeOz("protein_mappings.ids", m$universe)
  writeLines(m$id_open$ProteinId, "id_with_open.txt")
  writeLines(m$standard_annotation$ProteinId, "from_downloaded_db.txt")
  writeLines(m$universe$ProteinId, "universe.txt")
}
