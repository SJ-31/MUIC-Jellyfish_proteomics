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

# Gene Set enrichment analysis
gene_sets <- list(
  unknown_to_db = d$sample_tb %>%
    filter(inferred_by == "interpro" |
      inferred_by == "eggNOG" |
      grepl("[DT]", ProteinId)) %>%
    pluck("ProteinId"),
  has_mods = d$sample_tb %>%
    filter(ID_method == "open" | !is.na(mods)) %>%
    pluck("ProteinId")
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

# pwy <- group_pathways(d$sample_tb)
terms <- list(
  porin_activity = "GO:0015288",
  proteolysis = "GO:0006508",
  phospholipase_activity = "GO:0004620",
  pore_forming_activity = "GO:0140911",
  protein_folding = "GO:0006457",
  folding_chaperone = "GO:0140662",
  apoptotic_process = "GO:0006915",
  antioxidant_activity = "GO:0016209",
  metallopeptidase_activity = "GO:0008237",
  transporter_activity = "GO:0005215"
)
pwy <- group_go(d$sample_tb, wanted_terms = terms)

categories <- table(d$sample_tb$assigned_COG) %>%
  discard(., \(x) x < 100) %>%
  names()
category_lists <- group_list_from_tb(d$sample_tb, categories, "assigned_COG", "ProteinId")

grouping_col <- "GroupUP"
gene_sets <- c(gene_sets, pwy$grouped, category_lists)

by_intensity <- merge_lfq(d$sample_tb, "mean") %>%
  inner_join(., dplyr::select(d$sample_tb, c(ProteinId, {{ grouping_col }})))

fgsea_result_file <- glue("{fgsea_dir}/results.tsv")
fgsea_percolator_groups <- fgseaGroup(by_intensity,
  grouping_col = grouping_col,
  gene_sets = gene_sets
)

fgsea_sig <- fgsea_percolator_groups$fgsea$result %>% filter(padj < 0.05)

top <- fgsea_sig[head(order(NES, decreasing = TRUE), n = 15), pathway]

group_tb <- named_list2tb(fgsea_percolator_groups$groups)

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
  write_tsv(fgsea_sig, glue("{fgsea_dir}/formatted_results.tsv"))
  plots <- plotFgsea(
    fgsea_percolator_groups$groups,
    fgsea_percolator_groups$fgsea$ranked,
    fgsea_percolator_groups$fgsea$result
  )
  lapply(names(plots), \(x) {
    ggsave(glue("{fgsea_dir}/enrich_plot/{x}.png"), plots[[x]],
      create.dir = TRUE
    )
  })
} else {
  cat("", file = glue("{fgsea_dir}/no_results"))
}

# Note: The sizes are correctly the sizes of the ProteinGroups, not individual
# protein ids
fgsea_sig <- fgsea_percolator_groups$fgsea$result |> as_tibble()

go_id_map <- d$sample_tb |>
  separate_longer_delim(GO_IDs, ";") |>
  select(ProteinId, GO_IDs) |>
  filter(!is.na(GO_IDs))

filtered_pathways <- lapply(formatted$pathway, \(x) {
  current <- go_id_map |> filter(GO_IDs == x)
  name <- pathway_id_map[[x]]
  result <- list()
  result[[name]] <- d$sample_tb |> filter(ProteinId %in% current$ProteinId)
  result
}) |>
  unlist(recursive = FALSE)


for (n in names(filtered_pathways)) {
  nn <- str_replace_all(n, " ", "_")
  cur <- filtered_pathways[[n]] |>
    relocate(entry_name, pcoverage_align, .before = everything()) |>
    select(-c(ProteinGroupId, GroupSB))
  if (nrow(cur) > 0) {
    if (!dir.exists(glue("{fgsea_dir}/pathway_proteins"))) {
      dir.create(glue("{fgsea_dir}/pathway_proteins"))
    }
    write_tsv(cur, glue("{fgsea_dir}/pathway_proteins/{nn}.tsv"))
  }
}


TABLES$fgsea_sig_gt <- fgsea_sig |> gt()


save(c(TABLES, GRAPHS), fgsea_dir)
