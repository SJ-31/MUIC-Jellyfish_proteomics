if (!exists("SOURCED")) {
  source(paste0(getwd(), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
TABLES <- list()
GRAPHS <- list()
CACHE <- glue("{M$wd}/.cache")
library("enrichplot")

fgsea_dir <- glue("{M$outdir}/fgsea")
source(glue("{M$r_source}/KEGG_helpers.r"))

d <- get_go_data(M$data_w_cat_path,
  onto_path = M$ontologizer_path
)
perc_mods <- read_tsv(M$percolator_all_path) |>
  filter(ProteinId %in% d$sample_tb$ProteinId) |>
  group_by(ProteinId) |>
  summarise(mods = paste0(unique(mods |> discard(is.na)), collapse = ";")) |>
  filter(!mods == "")

d$sample_tb <- d$sample_tb |>
  left_join(perc_mods, by = join_by(ProteinId)) |>
  inner_join(M$lfq, by = join_by(ProteinId))

toxins <- d$sample_tb |>
  inner_join(read_tsv(M$toxin_map_path), by = join_by(GroupUP)) |>
  filter(!is.na(Group)) |>
  pluck("ProteinId")

# Gene Set enrichment analysis
gene_sets <- list(
  unknown_to_db = d$sample_tb %>%
    filter(inferred_by == "interpro" |
      inferred_by == "eggNOG" |
      grepl("[DT]", ProteinId)) %>%
    pluck("ProteinId"),
  has_mods = d$sample_tb %>%
    filter(ID_method == "open" | !is.na(mods)) %>%
    pluck("ProteinId"),
  toxins = toxins
)

kegg_metadata_file <- glue("{fgsea_dir}/kegg_metadata.tsv")
if (!file.exists(kegg_metadata_file)) {
  hh <- new.env()
  reticulate::source_python(glue("{args$python_source}/helpers.py"), envir = hh)
  kegg_entries <- fgsea_sig$pathway |> keep(\(x) str_detect(x, "[0-9]"))
  kegg_metadata <- hh$get_kegg_metadata(kegg_entries) |> as_tibble()
  write_tsv(kegg_metadata, kegg_metadata_file)
} else {
  kegg_metadata <- read_tsv(kegg_metadata_file)
}

pwy <- group_pathways(d$sample_tb)

categories <- table(d$sample_tb$assigned_COG) %>%
  discard(., \(x) x < 100) %>%
  names()
category_lists <- group_list_from_tb(d$sample_tb, categories, "assigned_COG", "ProteinId")

grouping_col <- "GroupUP"
gene_sets <- c(gene_sets, pwy$grouped, category_lists)
by_intensity <- merge_lfq(d$sample_tb, "mean") %>%
  inner_join(., dplyr::select(d$sample_tb, c(ProteinId, {{ grouping_col }})))

fgsea_result_file <- glue("{fgsea_dir}/results.tsv")
if (!file.exists(fgsea_result_file)) {
  fgsea_percolator_groups <- fgseaGroup(by_intensity,
    grouping_col = grouping_col,
    gene_sets = gene_sets
  )

  fgsea_sig <- fgsea_percolator_groups$fgsea$result %>% filter(padj < 0.05)

  keggid2pathway_name <- tb2named_list(kegg_metadata, "entry", "name")
  fgsea_sig <- fgsea_sig |>
    mutate(
      new_pathway =
        map_chr(pathway, \(x) ifelse(x %in% names(keggid2pathway_name), keggid2pathway_name[[x]], x))
    ) |>
    distinct(new_pathway, .keep_all = TRUE)

  renamed_groups <- local({
    g <- fgsea_percolator_groups$groups
    g <- g[names(g) %in% fgsea_sig$pathway]
    names(g) <- fgsea_sig$new_pathway
    g
  })
  fgsea_sig <- fgsea_sig |>
    mutate(pathway = new_pathway) |>
    select(-new_pathway)

  top <- fgsea_sig[head(order(NES, decreasing = TRUE), n = 15), pathway]

  group_tb <- named_list2tb(fgsea_percolator_groups$groups) |>
    mutate(name = map_chr(name, \(x) lget(keggid2pathway_name, x, x)))

  rank_tb <- fgsea_percolator_groups$fgsea$ranked |>
    as_tibble() |>
    mutate(group = names(fgsea_percolator_groups$fgsea$ranked), rank = rank(-value))

  formatted <- fgsea_sig |>
    as_tibble() |>
    inner_join(group_tb, by = join_by(x$pathway == y$name)) |>
    inner_join(rank_tb, by = join_by(x$value == y$group))

  top_enriched <- formatted |>
    filter(pathway %in% top) |>
    mutate(padj = -log(padj)) |>
    ggplot(aes(x = rank, y = NES, color = pathway, alpha = padj)) +
    geom_point(shape = 3, size = 4) +
    theme_bw() +
    M$default_theme +
    guides(
      color = guide_legend("Protein subset\n(In order of appearance from top)"),
      alpha = guide_legend("-Log adjusted p-value"),
    ) +
    ylab("Normalized Enrichment Score (NES)") +
    xlab("Intensity rank") +
    scale_x_continuous(
      breaks = seq(0, max(formatted$rank), by = 500)
    ) +
    scale_y_continuous(
      breaks = seq(4, max(formatted$NES), by = 0.2)
    ) +
    scale_color_discrete(breaks = top)
  GRAPHS$top_enriched <- top_enriched
  attr(GRAPHS$top_enriched, "width") <- 16


  if (!dir.exists(fgsea_dir)) {
    dir.create(fgsea_dir)
  }

  if (nrow(fgsea_percolator_groups$fgsea$result) != 0) {
    write_tsv(fgsea_percolator_groups$fgsea$result, glue("{fgsea_dir}/results.tsv"))
    write_tsv(fgsea_sig, glue("{fgsea_dir}/formatted_results"))
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
  # Note: The sizes are correctly the sizes of the ProteinGroups, not individual
  # protein ids
  fgsea_sig <- fgsea_percolator_groups$fgsea$result |> as_tibble()
} else {
  fgsea_sig <- read_tsv(fgsea_result_file)
}

kegg_sig <- fgsea_sig |> inner_join(kegg_metadata, by = join_by(x$pathway == y$entry))

TABLES$fgsea_sig_gt <- fgsea_sig |> gt()


save(c(TABLES, GRAPHS), fgsea_dir)
