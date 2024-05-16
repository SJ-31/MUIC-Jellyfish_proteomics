# Gene Set enrichment analysis

#' Ranked gene set enrichment analysis
#'
#' @description
#' Rank proteins in tb according to a specified column, in
#' descending fashion
#' @param quant the column to rank proteins on
fgseaWrapper <- function(quant, tb, gene_sets, id_col = "ProteinId") {
  ranked <- tb %>%
    dplyr::select(all_of(c(id_col, quant))) %>%
    dplyr::filter(!is.na(!!quant)) %>%
    column_to_rownames(var = id_col) %>%
    (\(x) {
      y <- as.vector(x[[quant]])
      names(y) <- rownames(x)
      return(y)
    }) %>%
    sort(., decreasing = TRUE)
  # Ranks are descending, so most abundant first
  fgsea <- fgsea::fgsea(pathways = gene_sets, ranked, scoreType = "pos")
  # Recommended to switch to "pos"
  return(list(result = fgsea, ranked = ranked))
}

# See here for more plotting
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
plotFgsea <- function(gene_sets, ranked_list, fgsea_result, p_cutoff = 0.05) {
  fgsea_plots <- list()
  for (n in seq_len(nrow(fgsea_result))) {
    if (fgsea_result[n, ]$padj < p_cutoff) {
      set <- fgsea_result[n, ]$pathway
      padjust <- round(fgsea_result[n, ]$padj, 6)
      fgsea_plots[[set]] <- plotEnrichment(gene_sets[[set]], ranked_list) +
        labs(title = glue("Set: {set}, adjusted p-value = {padjust}"))
    }
  }
  return(fgsea_plots)
}

fgseaGroup <- function(tb, grouping_col, gene_sets) {
  group_map <- hash::hash()
  group_map[tb$ProteinId] <- tb[[grouping_col]]
  grouped_intensity <- tb %>%
    group_by(!!as.symbol(grouping_col)) %>%
    nest() %>%
    mutate(
      mean_intensity = map_dbl(data, \(x) mean(x$log_intensity, na.rm = TRUE)),
      sd_intensity = map_dbl(data, \(x) sd(x$log_intensity, na.rm = TRUE)),
      members = lapply(data, \(x) x$ProteinId),
      size = map_dbl(data, \(x) nrow(x))
    ) %>%
    filter(!is.na(mean_intensity) & !is.na(!!as.symbol(grouping_col)))

  gene_sets_groups <- lapply(gene_sets, \(x) mapUnique(x, group_map))
  f <- fgseaWrapper("mean_intensity", grouped_intensity, gene_sets_groups, id_col = "Group")
  return(list(fgsea = f, groups = gene_sets_groups))
}

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

fgsea_percolator_groups <- fgseaGroup(by_intensity, grouping_col = "Group", gene_sets = gene_sets)
sig <- fgsea_percolator_groups$fgsea$result %>% filter(padj < 0.05)

plots <- plotFgsea(
  fgsea_percolator_groups$groups,
  fgsea_percolator_groups$fgsea$ranked,
  fgsea_percolator_groups$fgsea$result
)

write_tsv(fgsea_percolator_groups$fgsea$result, glue("{OUTDIR}/fgsea/results.tsv"))
lapply(names(plots), \(x) {
  ggsave(glue("{OUTDIR}/fgsea/{x}.png", plots[[x]]))
})
