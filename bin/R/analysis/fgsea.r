if (!exists("SOURCED")) {
  source(paste0(getwd(), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

source(glue("{M$r_source}/KEGG_helpers.r"))

d <- goData(M$combined_results,
  onto_path = M$ontologizer_path
)
# Gene Set enrichment analysis
gene_sets <- list(
  unknown_to_db = d$sample_tb %>%
    filter(inferred_by == "interpro" |
      inferred_by == "eggNOG" |
      grepl("[UDT]", ProteinId)) %>%
    pluck("ProteinId"),
  has_mods = d$sample_tb %>%
    filter(ID_method == "open" | !is.na(Mods)) %>%
    pluck("ProteinId"),
  toxins = names(getToxinProteins(d$prot_go_map$sample))
)
pwy <- groupPathways(d$sample_tb)
categories <- table(d$sample_tb$category) %>%
  discard(., \(x) x < 100) %>%
  names()
category_lists <- groupListFromTb(d$sample_tb, categories, "category", "ProteinId")

gene_sets <- c(gene_sets, pwy$grouped, category_lists)
by_intensity <- mergeLfq(d$sample_tb, "mean") %>%
  inner_join(., dplyr::select(d$sample_tb, c(ProteinId, Group)))

fgsea_percolator_groups <- fgseaGroup(by_intensity,
  grouping_col = "Group",
  gene_sets = gene_sets
)
fgsea_percolator_groups$fgsea$result <- fgsea_percolator_groups$fgsea$result %>% filter(padj < 0.05)


fgsea_dir <- glue("{M$outdir}/figures/fgsea")
if (!dir.exists(fgsea_dir)) {
  dir.create(fgsea_dir)
}
if (nrow(fgsea_percolator_groups$fgsea$result) != 0) {
  write_tsv(fgsea_percolator_groups$fgsea$result, glue("{fgsea_dir}/results.tsv"))
  plots <- plotFgsea(
    fgsea_percolator_groups$groups,
    fgsea_percolator_groups$fgsea$ranked,
    fgsea_percolator_groups$fgsea$result
  )
  lapply(names(plots), \(x) {
    ggsave(glue("{fgsea_dir}/{x}.png"), plots[[x]])
  })
} else {
  cat("", file = glue("{fgsea_dir}/no_results"))
}


# Visualize KEGG
GET_KEGG <- FALSE
if (GET_KEGG) {
  to_kegg <- select(M$data, -Group) |> inner_join(by_intensity, by = join_by(ProteinId))

  GENES <- flattenBy(to_kegg$KEGG_Genes, ";") |> unique()
  KO <- flattenBy(to_kegg$KEGG_ko, ";") |>
    unique() |>
    map_chr(\(x) glue("ko:{x}"))

  pathway_table <- flattenBy(to_kegg$KEGG_Pathway, ";") |>
    table() |>
    sort()
  module_table <- flattenBy(to_kegg$KEGG_Module, ";") |>
    table() |>
    sort()
  PATHWAYS <- names(pathway_table)
  MODULES <- names(module_table)

  CACHE <- glue("{M$wd}/.cache")


  pathway_graphs <- lapply(names(pathway_table), \(x) {
    pathway <- NULL
    try(pathway <- ggkegg::pathway(x, directory = CACHE))
    pathway
  }) |> `names<-`(names(pathway_table))

  pathway_tbs <- lapply(discard(pathway_graphs, is.null), \(x) as_tibble(tidygraph::activate(x, nodes)))

  pathway_data <- tibble(name = names(pathway_tbs)) |>
    mutate(
      completeness = map_dbl(pathway_tbs, pathway_completeness),
      title = map_chr(name, get_pathway_title),
      n_edges = map_dbl(names(pathway_tbs), \(x) length(E(pathway_graphs[[x]])))
    ) |>
    arrange(desc(completeness))

  rm(to_kegg)
}
