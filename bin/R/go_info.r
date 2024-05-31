if (sys.nframe() == 0) {
  library("optparse")
  library("glue")
  library("tidyverse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"))
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-o", "--outdir"))
  args <- parse_args(parser)
  source(glue("{args$r_source}/GO_helpers.r"))
  goVector(args$input, go_column = "GO_IDs", unique = TRUE) |>
    goInfoTb() |>
    write_tsv(glue("{outdir}/all_go_info.tsv"))
}
