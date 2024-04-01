library(tidyverse)

# Temporary fixes for bugs in `sort_interpro` script
# - Removes duplicate ids from `matchedPeptideIds` column
# - Removes `query` column
cleanDuplicateIds <- function(tb) {
  if (!purrr::pluck_exists(tb, "ProteinId")) {
    ids <- tb$query
  } else {
    ids <- tb$ProteinId
  }
  cleaned_ids <- purrr::map2_chr(
    ids, tb$matchedPeptideIds,
    \(x, y) {
      if (!is.na(y) && grepl(pattern = x, y)) {
        cleaned <- str_split_1(y, ";") %>%
          discard(., \(z) z == x) %>%
          paste0(., collapse = ";")
        return(cleaned)
      } else {
        return(y)
      }
    }
  )
  tb$matchedPeptideIds <- cleaned_ids
  return(tb)
}

fixSep <- function(tb) {
  fixed <- tb %>% mutate(GO = map_chr(GOs,
                                      \(x) str_replace_all(x, ",", ";")))
}


if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"))
  parser <- add_option(parser, c("-r", "--rename"), default = FALSE, action = "store_true")
  parser <- add_option(parser, c("-c", "--check_go"), default = FALSE, action = "store_true")
  args <- parse_args(parser)
  tb <- read_tsv(args$input)
  if (!args$rename) {
    ## if (purrr::pluck_exists(tb, "query")) {
    ##   tb <- dplyr::select(tb, -query)
    ## }
    write_tsv(cleanDuplicateIds(tb), args$input)
  } else {
    write_tsv(dplyr::rename(tb, inferred_by = Anno_method), args$input)
  }
}
