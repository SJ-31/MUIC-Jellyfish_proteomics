library("glue")
library("tidyverse")
library("optparse")

aggregate_embeddings <- function(args) {
  original <- read_tsv(args$combined_results)
  nrow_before <- nrow(original)
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    args$embeddings,
    args$distances
  )
  print("Embedding data received")
  py_clusters <- new.env()
  reticulate::source_python(glue("{args$python_source}/clustering.py"),
    envir = py_clusters
  )
  dist <- e$cosine
  clusters <- local({
    f <- hclust_sk(dist, 0.1, "average", labels_only = FALSE)
    linkage_matrix <- py_clusters$linkageMatrix(f$fitted)
    py_clusters$saveDendogram(linkage_matrix,
      glue("{args$outdir}/dendogram.png"),
      cutoff = 0.1
    )
    f$labels
  })
  print("Clusters obtained")
  e$metadata <- merge_clusters(clusters, e$metadata, "ProteinId")
  if (nrow_before != nrow(e$metadata)) {
    print("Rows changed")
    print(glue("nrow before: {nrow_before}"))
    print(glue("nrow now: {nrow(e$metadata)}"))
  }
  merged <- original |>
    filter(!ProteinId %in% e$metadata$ProteinId) |>
    bind_rows(e$metadata)
  # Fix is just to add the proteins not in the previous
  write_tsv(merged, args$combined_results)
  if (!is.null(args$ontologizer_path)) {
    aggregate_helper(e$metadata, "cluster", args$ontologizer_path, args$go_path, args$outdir)
  }
}

aggregate_helper <- function(tb, grouping_col, ontologizer_path, go_path, outdir) {
  nested <- aggregate_metadata(tb, grouping_col)
  slims <- nested$GO_slims %>%
    lapply(., ids_into_ontology) %>%
    purrr::reduce(., merge_lists)
  slim_tb <- nested %>%
    select(!!grouping_col, size) %>%
    ungroup() %>%
    mutate(
      !!grouping_col := as.vector(!!as.symbol(grouping_col)),
      slim_CC = slims$CC,
      slim_BP = slims$BP,
      slim_MF = slims$MF
    ) %>%
    mutate(across(contains("slim_"), \(x) dplyr::na_if(x, "")))
  write_tsv(slim_tb, glue("{outdir}/{grouping_col}-aggregated_slims.tsv"))
  nested <- enrich_groups(tb, nested, ontologizer_path, go_path, group_name = grouping_col)
  save_aggregated(nested, glue("{outdir}/{grouping_col}-aggregated.tsv"))
}


if (sys.nframe() == 0) {
  parser <- OptionParser()
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-c", "--combined_results"))
  parser <- add_option(parser, c("-s", "--sample_name"))
  parser <- add_option(parser, c("-o", "--ontologizer_path"), default = NULL)
  parser <- add_option(parser, c("-p", "--python_source"))
  parser <- add_option(parser, c("-g", "--go_path"))
  parser <- add_option(parser, c("-e", "--embeddings"))
  parser <- add_option(parser, "--outdir")
  parser <- add_option(parser, c("-d", "--distances"))
  parser <- add_option(parser, "--aggregate", action = "store_true", default = FALSE)
  args <- parse_args(parser)
  source(glue("{args$r_source}/GO_helpers.r"))
  source(glue("{args$r_source}/helpers.r"))
  source(glue("{args$r_source}/analysis/metric_functions.r"))
  source(glue("{args$r_source}/cluster_helpers.r"))
  source(glue("{args$r_source}/rrvgo_modified.r"))
  source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
  tryCatch(
    expr = aggregate_embeddings(args),
    error = \(cnd) {
      last_error <- reticulate::py_last_error()
      message("Python error: ", last_error$type, "\n", last_error$value, "\n", last_error$traceback)
    }
  )
  if (args$aggregate) {
    aggregate_helper(
      read_tsv(args$combined_results), "Group",
      args$ontologizer_path, args$go_path, args$outdir
    )
    aggregate_helper(
      read_tsv(args$combined_results), "GroupUP",
      args$ontologizer_path, args$go_path, args$outdir
    )
  }
}
