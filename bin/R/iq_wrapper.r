library(tidyverse)
library(glue)
library(iq)

main <- function(dlfq_input) {
  dlfq <- read_tsv(dlfq_input)
  sample_cols <- colnames(dlfq) %>%
    purrr::discard(grepl("protein|ion", .))
  names(sample_cols) <- paste0("maxlfq-", sample_cols)
  dlfq <- dlfq %>% separate_longer_delim(., "protein", ";")
  dlfq_sep <- map(sample_cols, \(x) {
    dplyr::select(dlfq, c("protein", "ion", x)) %>%
      dplyr::filter(!!as.symbol(x) > 0) %>%
      dplyr::mutate(sample = x) %>%
      dplyr::rename(intensity = !!x)
  }) %>% bind_rows()

  preprocessed <- iq::preprocess(as.data.frame(dlfq_sep),
    primary_id = "protein", secondary_id = "ion",
    sample_id = "sample",
    intensity_col = "intensity"
  )

  result <- iq::fast_MaxLFQ(preprocessed)$estimate %>%
    as.data.frame() %>%
    rownames_to_column(var = "ProteinId") %>%
    as_tibble() %>%
    rename(., all_of(sample_cols))
  return(result)
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"), type = "character", help = "input file in form of directlfq generic input")
  parser <- add_option(parser, c("-r", "--r_source"), type = "character", help = "path to r source directory")
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  iq_result <- main(args$input)
  write_tsv(iq_result, args$output)
}
