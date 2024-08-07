TABLES <- list()
GRAPHS <- list()
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


# ----------------------------------------
# Relation between identification method and protein categories
category_col <- "GO_category_CC"
# Will choose cellular location, because this has the largest impact on acquisition
# among GO terms
cur <- M$data |> select(ID_method, {{ category_col }})
id_methods <- unique(cur$ID_method)
result <- chisqNME(
  cur, id_methods,
  category_col, "ID_method", category_col, "ID_method",
  binary = FALSE
)

TABLES$id_method_cc_chi <- result$gt$chi
TABLES$id_method_cc_contingency <- result$gt$contingency

category_col <- "GO_category_MF"
cur <- inner_join(M$data, M$taxa_tb, by = join_by(ProteinId)) |> select(Class, {{ category_col }})
classes <- unique(cur$Class) |> discard(is.na)
result <- chisqNME(cur, classes, category_col, "Class", category_col, "Class",
  binary = FALSE
)

TABLES$class_mf_chi <- result$gt$chi
TABLES$class_mf_contingency <- result$gt$contingency

# ----------------------------------------
# Evaluate the denovo proteins belonging to chironex that
grouped <- M$data |> group_by(GroupUP)
has_denovo <- grouped |>
  filter(!is.na(MatchedPeptideIds) & grepl("D", MatchedPeptideIds)) |>
  inner_join(M$taxa_tb, by = join_by(ProteinId))

# with_taxa <- M$repr |> inner_join(M$taxa_tb, by = join_by(ProteinId)) |>

is_chironex <- has_denovo |> filter(Genus == "Chironex")

not_chironex <- has_denovo |> filter(Genus != "Chironex")
not_chironex$Genus |> table()



save(c(TABLES, GRAPHS), glue("{M$outdir}/subgroup_comparisons"))
