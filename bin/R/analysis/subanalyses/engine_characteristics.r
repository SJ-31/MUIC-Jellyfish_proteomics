if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
PALETTE <- "ggthemes::colorblind"
PALETTE2 <- "ggthemes::Classic_10_Medium"
library("ggplot2")
library("ggVennDiagram")
TABLES <- list()
GRAPHS <- list()

# File and path setup
open_search_engines <- c("metamorpheusGTPMD", "msfraggerGPTMD", "msfraggerGlyco")
percolator_all <- read_tsv(M$percolator_all)
ENGINES <- percolator_all$engine |> unique()
standard_search_engines <- ENGINES[!ENGINES %in% open_search_engines]


data <- M$data

ta <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = ta)
per_protein_alignment_metrics_file <- glue("{M$outdir}/per_protein_alignment_metrics.tsv")
if (file.exists(per_protein_alignment_metrics_file)) {
  per_protein_alignment_metrics <- read_tsv(per_protein_alignment_metrics_file)
} else {
  tracer <- ta$AlignmentTracer(M$aligned_peptides_path, M$peptide_map_path)
  per_protein_alignment_metrics <- tracer$run() |> as_tibble()
  write_tsv(per_protein_alignment_metrics, per_protein_alignment_metrics_file)
}

engine_alignment_metrics <- per_protein_alignment_metrics |> select(
  ProteinId,
  contains(ENGINES), -matches("count|unmatched")
)

GRAPHS$engine_peptide_coverage <- engine_alignment_metrics %>%
  pivot_longer(., cols = -ProteinId) %>%
  mutate(
    value = log(value * 100),
    name = map_chr(name, \(x) str_replace(x, "_coverage", "")),
    type = case_when(
      name %in% standard_search_engines ~ "standard",
      name %in% open_search_engines ~ "open search",
    )
  ) |>
  ggplot(aes(y = value, x = name, fill = name, color = type)) +
  geom_boxplot() +
  labs(title = "Density of per-engine peptide coverage") +
  ylab("log coverage (%)") +
  xlab("Engine name") +
  guides(fill = "none") +
  scale_fill_paletteer_d(PALETTE2) +
  scale_color_manual(values = c("standard" = "#4c4f69", "open search" = "#d20f39")) +
  scale_x_discrete(
    limits = c(standard_search_engines, open_search_engines)
  )

# Identify which engine, if any, is the best-performing
cov_list <- engine_alignment_metrics %>%
  select(-ProteinId) |>
  as.list()
names(cov_list) <- names(cov_list) |> map_chr(\(x) str_replace(x, "_coverage", ""))
ks <- kruskal.test(cov_list)
combos <- combn(names(cov_list), 2)
test_tb <- lapply(
  seq_len(ncol(combos)),
  \(x) {
    greater <- wilcox.test(cov_list[[combos[1, x]]],
      cov_list[[combos[2, x]]],
      alternative = "greater"
    )
    greater$data.name <- glue("{combos[1, x]} x {combos[2, x]}")
    greater$alternative <- glue("{combos[1, x]} greater")
    two_sided <- wilcox.test(
      cov_list[[combos[1, x]]],
      cov_list[[combos[2, x]]]
    )
    two_sided$data.name <- glue("{combos[1, x]} x {combos[2, x]}")
    two_sided$alternative <- glue("two sided")
    bind_rows(htest2tb(greater), htest2tb(two_sided))
  }
) %>% bind_rows()
test_tb <- test_tb %>%
  dplyr::select(-c(method, null)) %>%
  mutate(
    p_adjust = p.adjust(p_value),
    significant = map_dbl(p_adjust, \(x) x < 0.05)
  ) %>%
  rename(pair = data)
TABLES$engine_coverage_pairwise <- gt(test_tb)
TABLES$engine_coverage_pairwise_sig <- gt(conclude_one_sided(test_tb))

num_peptides_matched <- ta$get_engine_counts(M$percolator_all, data) |> as_tibble()
# #' Group engines that identify the same peptide groups using Jaccard distance
enginesXProtein <- num_peptides_matched %>%
  # select(-all_of(open_search_engines)) %>%
  tb_transpose()

engine_dist <- vegan::vegdist(enginesXProtein, method = "jaccard")
dist_longer <- engine_dist |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = "engine") |>
  pivot_longer(cols = -engine)
pcoa <- vegan::wcmdscale(engine_dist, eig = TRUE)$points %>%
  as.data.frame() |>
  rownames_to_column(var = "engine")
GRAPHS$engine_sim_jaccard <- dist_longer |> ggplot(aes(x = engine, y = name, fill = value)) +
  geom_tile() +
  ylab("Engine") +
  xlab("Engine") +
  scale_fill_paletteer_c("ggthemes::Orange-Gold", name = "Jaccard distance") +
  theme(axis.text.x = element_text(angle = 90))
GRAPHS$engine_sim_jaccard_biplot <- pcoa |> ggplot(aes(x = Dim1, y = Dim2, color = engine)) +
  geom_point()



standard_engine_tbs <- lapply(standard_search_engines, \(x) {
  percolator_all |> filter(engine == x)
}) %>%
  `names<-`(standard_search_engines)

#' Venn diagram for overlap
id_list <- lapply(standard_engine_tbs, \(x) x$ProteinId) %>%
  `names<-`(names(standard_engine_tbs))
venn <- ggVennDiagram(id_list, label = "none") +
  scale_fill_paletteer_c("ggthemes::Classic Blue")
GRAPHS$engine_venn <- venn

# ----------------------------------------
# Evaluate the contribution of each alignment type onto the protein
no_denovo <- read_tsv(glue("{M$wd}/results/ND_C_indra/{M$chosen_pass}/ND_C_indra_all_wcoverage.tsv"))
type_alignment_metrics <- per_protein_alignment_metrics |>
  select(
    ProteinId, contains(c("transcriptome", "unmatched", "denovo", "database")),
    -contains("count")
  ) |>
  inner_join(select(no_denovo, ProteinId, pcoverage_nmatch)) |>
  rename(database_nd = pcoverage_nmatch) |>
  inner_join(select(data, ProteinId, pcoverage_align)) |>
  rename(total = pcoverage_align)

type_alignment_metrics_long <- type_alignment_metrics |>
  pivot_longer(cols = -ProteinId) |>
  mutate(name = map_chr(\(x) str_replace("_coverage", "")))


save(c(GRAPHS, TABLES), glue("{M$outdir}/figures/engine_characteristics"))
