if (sys.nframe() == 0) {
  library("optparse")
  library("glue")
  library("tidyverse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"))
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-f", "--go_info_path"))
  parser <- add_option(parser, c("-p", "--python_source"))
  parser <- add_option(parser, c("-d", "--predefined_groups"))
  parser <- add_option(parser, c("-g", "--go_path"))
  parser <- add_option(parser, c("-n", "--n_groups"))
  parser <- add_option(parser, c("-j", "--parent_mapping"))
  parser <- add_option(parser, c("-o", "--outdir"))
  parser <- add_option(parser, c("-m", "--mode"))
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  source(glue("{args$r_source}/GO_helpers.r"))
  if (args$mode == "info") {
    d <- get_go_data(args$input)
    get_go_vec(d$sample_tb, go_column = "GO_IDs", unique = TRUE) |>
      go_info_tb() |>
      write_tsv(glue("{args$outdir}/all_go_info.tsv"))
  } else if (args$mode == "get_parents") {
    reticulate::source_python(glue("{args$python_source}/go_subset.py"))
    to_json(
      go_path = args$go_path, sample_path = args$input,
      go_info_path = args$go_info_path, outdir = args$outdir,
      predefined = args$predefined_groups,
      n_groups = as.integer(args$n_groups)
    )
  } else if (args$mode == "write_results") {
    reticulate::source_python(glue("{args$python_source}/go_subset.py"))
    result <- reticulate_show_error(find_go_parents(
      args$input, args$go_info_path, args$parent_mapping, args$predefined_groups
    )) |> lapply(as_tibble)
    new_name <- str_replace(args$input, "wcoverage", "wcategory")
    result[[1]] %>%
      dplyr::relocate(., contains("GO_category"), .after = "GO_slims") %>%
      write_tsv(., new_name)
    result[[2]] %>% write_tsv(., glue("{args$outdir}/missing_counts.tsv"))
  }
}
