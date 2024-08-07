library("optparse")
library("glue")

main <- function(args) {
  library(reticulate)
  source(glue("{args$r_source}/GO_helpers.r"))
  source(glue("{args$r_source}/DR_helpers.r"))
  source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
  if (!is.null(args$compare) && args$compare) {
    e <- embeddingData(args$combined_results,
      args$sample_name,
      args$embedding_path,
      args$dist_path,
      comparison_meta = args$uniprot_data
    )
  } else {
    e <- embeddingData(
      args$combined_results,
      args$sample_name,
      args$embedding_path,
      args$dist_path
    )
  }

  prefix <- ifelse(is.null(args$results_prefix), args$sample_name,
    args$results_prefix
  )
  result <- dr_wrapper(
    e, "ProteinId",
    args$figure_path,
    prefix,
    args$technique
  )
  title_str <- ifelse(!is.null(args$compare),
    glue("{args$sample_name} comparison"),
    glue("{args$sample_name} sample")
  )
  for (color in e$color) {
    label <- label_gen(args$technique, title_str, "")
    plot_dr(
      to_plot = result$to_plot,
      color_col = color,
      path = args$figure_path,
      technique = args$technique,
      labels = label,
      twod = TRUE,
      unwanted = c("unknown", "other")
    )
  }
}

if (sys.nframe() == 0 && length(commandArgs(TRUE))) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-f", "--figure_path"))
  parser <- add_option(parser, c("-s", "--sample_name"))
  parser <- add_option(parser, c("-u", "--uniprot_data"))
  parser <- add_option(parser, c("-t", "--technique"))
  parser <- add_option(parser, c("-e", "--embedding_path"))
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-d", "--dist_path"))
  parser <- add_option(parser, c("-p", "--python_source"))
  parser <- add_option(parser, "--results_prefix")
  parser <- add_option(parser, "--compare", action = "store_true")
  parser <- add_option(parser, c("-c", "--combined_results"))
  args <- parse_args(parser)
  main(args)
}
