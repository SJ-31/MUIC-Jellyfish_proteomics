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
alignment_types <- c("denovo", "transcriptome", "database", "unmatched_peptide")
standard_search_engines <- ENGINES[!ENGINES %in% open_search_engines]


data <- M$data
ta <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = ta)

# Get general alignment metrics
per_protein_alignment_metrics_file <- glue("{M$outdir}/per_protein_alignment_metrics.tsv")
if (file.exists(per_protein_alignment_metrics_file)) {
  per_protein_alignment_metrics <- read_tsv(per_protein_alignment_metrics_file)
} else {
  tracer <- ta$AlignmentTracer(M$aligned_peptides_path, M$peptide_map_path)
  per_protein_alignment_metrics <- tracer$run() |> as_tibble()
  write_tsv(per_protein_alignment_metrics, per_protein_alignment_metrics_file)
}

# Get file for evaluating cost of removing a type of alignment from data
per_protein_alignment_differences_file <- glue("{M$outdir}/per_protein_alignment_differences.tsv")
if (file.exists(per_protein_alignment_differences_file)) {
  per_protein_alignment_differences <- read_tsv(per_protein_alignment_differences_file)
} else {
  tracer <- ta$AlignmentTracer(M$aligned_peptides_path, M$peptide_map_path)
  per_protein_alignment_differences <- tracer$run(mode = "differences") |> as_tibble()
  write_tsv(per_protein_alignment_differences, per_protein_alignment_differences_file)
}


# ----------------------------------------
# Engine alignments
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

engine_alignment_metrics

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
TABLES$engine_coverage_pairwise_sig <- conclude_one_sided(test_tb) |> pairwise_conclusion2gt()

num_peptides_matched <- ta$get_engine_counts(M$percolator_all, data) |>
  as_tibble() |>
  distinct()

# #' Group engines that identify the same peptide groups using Jaccard distance

enginesXProtein <- num_peptides_matched %>%
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
# Evaluate the contribution of each engine onto the protein

engine_tb <- per_protein_alignment_differences |>
  select(ProteinId, contains(ENGINES)) |>
  select(-unmatched_peptide)
engine_longer <- engine_tb |>
  pivot_longer(cols = -ProteinId)

cov_list <- engine_tb |>
  select(-ProteinId) |>
  as.list()
cov_list$total <- per_protein_alignment_differences$total_coverage


cov_tests <- test_all_pairs(cov_list, wilcox.test, two_sided = TRUE) |>
  bind_rows(
    test_all_pairs(cov_list, \(x, y) wilcox.test(x, y, alternative = "less"), alternative_suffix = "less")
  ) |>
  mutate(p_adjust = p.adjust(p_value), significant = ifelse(p_adjust < 0.05, 1, 0))

# "x less than y" means that removing engine x's peptides has a greater impact on coverage compared to y
cov_test_conclusion <- conclude_one_sided(cov_tests)
if (nrow(cov_test_conclusion) > 0) {
  cov_test_conclusion <- filter(!is.na(conclusion))
  TABLES$coverage_impact_conclusion <- gt(cov_test_conclusion)
} else {
  TABLES$no_significant_difference_when_removing_engine_peptides <- 0
}


GRAPHS$engines_removed <- engine_longer |> ggplot(aes(y = value, fill = name)) +
  geom_boxplot() +
  scale_fill_paletteer_d(PALETTE2) +
  guides(fill = guide_legend(title = "Engine")) +
  ylab("Coverage (%) if engine's peptides were removed")

# ----------------------------------------
# Comparison of engine peptide characteristics
pepmap <- read_tsv(M$peptide_map_path)
all_mapped_scans <- read_tsv(M$mapped_scan_path)

joined <- inner_join(pepmap, all_mapped_scans,
  by = join_by(x$peptideIds == y$base_peptide, engine)
)

# Sets of peptides per engine
peptide_sets <- lapply(unique(joined$engine), \(x) {
  pepmap |>
    filter(engine == x)
}) |>
  `names<-`(unique(joined$engine))

test_wrapper <- function(var) {
  sets <- col_from_tb_list(peptide_sets, var) |> discard(\(x) all(is.na(x)))
  kruskal <- kruskal.test(sets)
  if (kruskal$p.value > 0.05) {
    warning(glue("Pairwise tests for peptide variable `{var}` not run, kruskal test not significant"))
    return(htest2tb(kruskal))
  }
  tests <- test_all_pairs(sets, \(x, y) wilcox.test(x, y, alternative = "greater"), "greater") |>
    bind_rows(test_all_pairs(sets, wilcox.test, two_sided = TRUE)) |>
    mutate(
      p_adjust = p.adjust(p_value),
      significant = ifelse(p_adjust < 0.05, 1, 0)
    )
  conclude_one_sided(tests)
}

# Compare mass, length and ???
masses <- col_from_tb_list(peptide_sets, "mass") |> discard(\(x) all(is.na(x)))
lengths <- col_from_tb_list(peptide_sets, "length") |> discard(\(x) all(is.na(x)))

mass_plot <- gg_numeric_dist(lapply(masses, \(x) log(x)), method = "boxplot") +
  ylab("log mw") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  guides(color = "none")

length_plot <- gg_numeric_dist(lapply(lengths, \(x) log(x)), method = "boxplot") +
  ylab("Peptide length") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  guides(color = "none")


TABLES$mass_test_result <- test_wrapper("mass") |> pairwise_conclusion2gt()
TABLES$length_test_result <- test_wrapper("length") |> pairwise_conclusion2gt()


# ----------------------------------------
#' Engine peptide characteristics, per-protein, paired tests

matched_peptide_map <- {
  temp <- M$data |>
    filter(!is.na(MatchedPeptideIds)) |>
    select(ProteinId, MatchedPeptideIds) |>
    separate_longer_delim(MatchedPeptideIds, ";")
  temp
}


# Map any matched UP peptides back into DBPs
pepmap_joined <- left_join(pepmap, matched_peptide_map,
  by = join_by(x$ProteinId == y$MatchedPeptideIds),
  relationship = "many-to-many"
) |>
  filter(!is.na(ProteinId.y) | grepl("P", ProteinId)) |>
  mutate(
    ProteinId = map2_chr(ProteinId, ProteinId.y, \(x, y) {
      ifelse(str_detect("P", x) && is.na(y), x, y)
    })
  ) |>
  select(-ProteinId.y) |>
  filter(!is.na(engine))

#' Tests engine peptide metrics on a per-protein basis,
#' rather than by engine groups
#' @param summarize_fn A function that aggregates the peptide variable `variable`
#' into a single value for each peptide set per engine
compare_engine_peptide_metrics <- function(
    peptide_map,
    variable = "mass",
    summarize_fn = \(x) summarize(x, mass = mean(mass))) {
  combos <- combn(unique(peptide_map$engine), 2)

  compare_pair <- function(pair) {
    filtered <- filter(peptide_map, engine %in% pair)
    cur_engines <- filtered |>
      group_by(ProteinId) |>
      filter(length(unique(engine)) == 2) |>
      group_by(ProteinId, engine) |>
      summarize_fn() |>
      ungroup() |>
      pivot_wider(names_from = engine, values_from = !!as.symbol(variable))
    x <- cur_engines[[pair[1]]]
    y <- cur_engines[[pair[2]]]
    p_string <- glue("{pair[1]} x {pair[2]}")
    t1 <- wilcox.test(x, y, paired = TRUE) |> htest2tb(data.name = p_string)
    t2 <- wilcox.test(x, y,
      paired = TRUE,
      alternative = "greater"
    ) |>
      htest2tb(data.name = p_string, alternative = glue("{pair[1]} greater"))
    bind_rows(t1, t2)
  }
  lapply(seq_len(ncol(combos)), \(i) compare_pair(combos[, i])) |>
    bind_rows() |>
    rename(pair = data)
}


mass_tests <- compare_engine_peptide_metrics(pepmap_joined) |>
  get_adjusted_p()
mass_test_conclusions <- conclude_one_sided(mass_tests)
length_tests <- compare_engine_peptide_metrics(pepmap_joined,
  "length",
  summarize_fn = \(x) summarize(x, length = mean(length))
) |>
  get_adjusted_p()
length_test_conclusions <- conclude_one_sided(length_tests)


# ----------------------------------------
#' Engine peptide overlaps
peps_list <- col_from_tb_list(peptide_sets, "peptideIds") |> discard(\(x) all(is.na(x)))
standard_peps <- peps_list[!names(peps_list) %in% open_search_engines]


overlap_results <- ta$get_peptide_overlap(M$peptide_map_path) |> `names<-`(c("overlap_df", "subsets"))
overlap_tb <- overlap_results$overlap_df |> as_tibble()
GRAPHS$engine_sim_overlap <- overlap_tb |> ggplot(aes(x = first, y = second, fill = overlap)) +
  geom_tile() +
  geom_text(aes(label = round(overlap, 2))) +
  ylab("Engine") +
  xlab("Engine") +
  scale_fill_paletteer_c("ggthemes::Classic Area Green", name = "Overlap coefficient") +
  theme(axis.text.x = element_text(angle = 90))

if (!dir.exists(glue("{M$outdir}/engine_characteristics"))) {
  dir.create(glue("{M$outdir}/engine_characteristics"))
}
if (length(overlap_results$subsets) == 0) {
  cat("no subsets",
    file = glue("{M$outdir}/engine_characteristics/peptide_subsets.txt")
  )
} else {
  cat(overlap_results$subsets,
    file = glue("{M$outdir}/engine_characteristics/peptide_subsets.txt")
  )
}
overlap_tb <- overlap_tb |> distinct(smaller, overlap, .keep_all = TRUE)
TABLES$overlap_coefficient <- overlap_tb


save(c(GRAPHS, TABLES), glue("{M$outdir}/engine_characteristics"))
