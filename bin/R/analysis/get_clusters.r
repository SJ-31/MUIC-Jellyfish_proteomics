source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/rrvgo_modified.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))

e <- embeddingData(
  args$combined_results,
  args$sample_name,
  PROTTRANS_EMBD,
  PROTTRANS_DIST
)

if (!file.exists(glue("{OUTDIR}/aggregated_clusters.tsv"))) {
  py_clusters <- new.env()
  reticulate::source_python(glue("{args$python_source}/clustering.py"),
    envir = py_clusters
  )
  dist <- e$cosine
  clusters <- local({
    f <- hclustSk(dist, 0.1, "average", labels_only = FALSE)
    linkage_matrix <- py_clusters$linkageMatrix(f$fitted)
    py_clusters$saveDendogram(linkage_matrix,
      glue("{OUTDIR}/dendogram.png"),
      cutoff = 0.1
    )
    f$labels
  })
  e$metadata <- mergeClusters(clusters, e$metadata, "ProteinId")
  nested <- aggregateMetadata(e$metadata, "cluster")
  slims <- nested$GO_slims %>%
    lapply(., idsIntoOntology) %>%
    purrr::reduce(., mergeLists)
  slim_tb <- nested %>%
    select(cluster, size) %>%
    ungroup() %>%
    mutate(
      cluster = as.vector(cluster),
      slim_CC = slims$CC,
      slim_BP = slims$BP,
      slim_MF = slims$MF
    ) %>%
    mutate(across(contains("slim_"), \(x) dplyr::na_if(x, "")))
  nested <- enrichGroups(
    e$metadata, nested,
    glue("{tools}/Ontologizer.jar"), args$go_path
  )
  saveAggregated(nested, glue("{OUTDIR}/aggregated_clusters.tsv"))
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
