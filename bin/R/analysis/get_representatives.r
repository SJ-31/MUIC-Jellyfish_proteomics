library("tidyverse")
#' Finds a representative protein from Protein Groups
#'
#' @description
#' "Protein Group" is defined here as proteins that all share the same set of peptides
#' The representative is the member with the highest *coverage in the set, unless
#' it is found to be a protein fragment (determined from the header), or it is the
#' the shortest protein in the set
get_representatives <- function(tb) {
  if (nrow(tb) == 1) {
    return(tb)
  }
  complete <- tb |> filter(!grepl("fragment|partial", header, ignore.case = TRUE))
  if (nrow(complete) == 0) {
    to_find <- tb
  } else {
    to_find <- complete
  }
  return(arrange(to_find, desc(pcoverage_align)) |> slice(1))
}


if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"), type = "character", action = "store")
  parser <- add_option(parser, c("-o", "--output"), type = "character", action = "store")
  args <- parse_args(parser)
  nested <- read_tsv(args$input) |>
    group_by(GroupUP) |>
    nest() |>
    mutate(size = map_dbl(data, nrow)) |>
    arrange(desc(size))
  lapply(nested$data, get_representatives) |>
    bind_rows() |>
    write_tsv(args$output)
}
