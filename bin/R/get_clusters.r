library("glue")
library("tidyverse")
library("optparse")

aggregateEmbeddings <- function(args) {
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
    f <- hclustSk(dist, 0.1, "average", labels_only = FALSE)
    linkage_matrix <- py_clusters$linkageMatrix(f$fitted)
    py_clusters$saveDendogram(linkage_matrix,
      glue("{args$outdir}/dendogram.png"),
      cutoff = 0.1
    )
    f$labels
  })
  print("Clusters obtained")
  e$metadata <- mergeClusters(clusters, e$metadata, "ProteinId")
  write_tsv(e$metadata, args$combined_results)
  aggregateHelper(e$metadata, "cluster", args$ontologizer_path, args$go_path, args$outdir)
  ## nested <- aggregateMetadata(e$metadata, "cluster")
  ## slims <- nested$GO_slims %>%
  ##   lapply(., idsIntoOntology) %>%
  ##   purrr::reduce(., mergeLists)
  ## slim_tb <- nested %>%
  ##   select(cluster, size) %>%
  ##   ungroup() %>%
  ##   mutate(
  ##     cluster = as.vector(cluster),
  ##     slim_CC = slims$CC,
  ##     slim_BP = slims$BP,
  ##     slim_MF = slims$MF
  ##   ) %>%
  ##   mutate(across(contains("slim_"), \(x) dplyr::na_if(x, "")))
  ## nested <- enrichGroups(
  ##   e$metadata, nested,
  ##   glue("{tools}/Ontologizer.jar"), args$go_path
  ## )
  ## saveAggregated(nested, glue("{OUTDIR}/aggregated_clusters.tsv"))
}

aggregateHelper <- function(tb, grouping_col, ontologizer_path, go_path, outdir) {
  nested <- aggregateMetadata(tb, grouping_col)
  slims <- nested$GO_slims %>%
    lapply(., idsIntoOntology) %>%
    purrr::reduce(., mergeLists)
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
  nested <- enrichGroups(tb, nested, ontologizer_path, go_path, group_name = grouping_col)
  saveAggregated(nested, glue("{outdir}/{grouping_col}-aggregated.tsv"))
}


if (sys.nframe() == 0) {
  parser <- OptionParser()
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-c", "--combined_results"))
  parser <- add_option(parser, c("-s", "--sample_name"))
  parser <- add_option(parser, c("-o", "--ontologizer_path"))
  parser <- add_option(parser, c("-p", "--python_source"))
  parser <- add_option(parser, c("-g", "--go_path"))
  parser <- add_option(parser, c("-e", "--embeddings"))
  parser <- add_option(parser, "--outdir")
  parser <- add_option(parser, c("-d", "--distances"))
  args <- parse_args(parser)
  source(glue("{args$r_source}/GO_helpers.r"))
  source(glue("{args$r_source}/helpers.r"))
  source(glue("{args$r_source}/analysis/metric_functions.r"))
  source(glue("{args$r_source}/cluster_helpers.r"))
  source(glue("{args$r_source}/rrvgo_modified.r"))
  source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
  tryCatch(
    expr = aggregateEmbeddings(args),
    error = \(cnd) {
      last_error <- reticulate::py_last_error()
      message("Python error: ", last_error$type, "\n", last_error$value, "\n", last_error$traceback)
    }
  )
  aggregateHelper(
    read_tsv(args$combined_results), "Group",
    args$ontologizer_path, args$go_path, args$outdir
  )
}

# TODO: Parse the lineages better
# sample <- e$metadata$lineage %>%
#   discard(is.na) %>%
#   index(1, 5)
# sample
# # summariseLineages <- function(lineage_vector) {
# longest <- 0
# taxa_splits <- lapply(sample, \(x) {
#   splits <- str_split_1(x, ";")
#   if (length(splits) > longest) {
#     longest <<- length(splits)
#   }
#   splits
# })
# m <- matrix(NA, ncol = longest, nrow = length(taxa_splits))
# for (n in seq_along(taxa_splits)) {
#   m[n, ] <- c(taxa_splits[[n]], rep(NA, longest - length(taxa_splits[[n]])))
# }
