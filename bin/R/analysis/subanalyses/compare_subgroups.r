#' Are there any meaningful differences between the sub-groups of proteins
#' identified in the data?
#'


#' Try to predict what an unannotated protein does based on what the members
#' of its group does (e.g. its GO_IDs and KEGG)

data <- M$data


cluster_aggregated <- read_tsv(glue("{M$path}/Analysis/Aggregated/cluster-aggregated.tsv"))
singleton_clusters <- cluster_aggregated |> filter(size == 1)
# singleton_groups <-
unknown <- data |> filter((is.na(GO_IDs) | is.na(header)))

unknown_single_by_cluster <- unknown |> filter(cluster %in% singleton_clusters$cluster)
unknown_others <- unknown |> filter(!ProteinId %in% unknown_single_by_cluster$ProteinId)

toxins <- cluster_aggregated |>
  filter(cluster %in% unknown_others$cluster) |>
  filter(if_any(contains("GO"), \(x) grepl("toxin", x)))

unknown_others |>
  filter(cluster %in% toxins$cluster) |>
  pluck("header")
