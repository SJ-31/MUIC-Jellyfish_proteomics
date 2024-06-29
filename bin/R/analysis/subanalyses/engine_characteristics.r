if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
PALETTE <- "ggthemes::colorblind"
library("ggplot2")
library("ggVennDiagram")
TABLES <- list()
GRAPHS <- list()

get_all_intervals <- function(ids) {
  get_id_intervals <- function(id) {
    lapply(ENGINES, \(x) py$tracer$engine_alignment_intervals(id, x)) %>%
      `names<-`(ENGINES)
  }
  lapply(ids, get_id_intervals) %>% `names<-`(ids)
}

interval_coverage <- function(interval_list, length, id) {
  sum_interval <- function(intervals) {
    if (length(intervals) == 0) {
      return(0)
    }
    s <- sum(map_dbl(intervals, \(x) x[2] - x[1]))
    s / length
  }
  lmap(
    interval_list,
    \(x) tibble(engine = names(x), ProteinId = id, coverage = sum_interval(x[[1]]))
  ) %>%
    bind_rows() %>%
    pivot_wider(names_from = engine, values_from = coverage)
}

# File and path setup
pass <- "1-First_pass"
PERCOLATOR_DIR <- c(glue("{M$path}/{pass}/Percolator"), glue("{M$path}/{pass}/Open_search/Percolator"))
open_search_engines <- c("metamorpheusGTPMD", "msfraggerGPTMD", "msfraggerGlyco")
percolator_files <- list.files(PERCOLATOR_DIR,
  pattern = "_percolator_proteins.tsv", full.names = TRUE
)
alignments <- alignmentData(M$path, "first")
combined_results <- M$data


percolator_all <- read_tsv(M$percolator_all)
ENGINES <- unique(percolator_all$engine)
percolator_tibbles <- list()
# Engine distribution analysis
for (e in ENGINES) {
  percolator_tibbles[[e]] <- percolator_all |> filter(engine == e)
}

standard_search <- combined_results %>% filter(ID_method == "standard")

normal_engine_tbs <- percolator_tibbles[!names(percolator_tibbles) %in% open_search_engines]


num_peptides_matched <- apply(combined_results, 1, \(x) {
  protein_id <- x["ProteinId"]
  ids_to_match <- protein_id
  matched_peptides <- x["MatchedPeptideIds"]
  if (!is.na(matched_peptides)) {
    ids_to_match <- c(ids_to_match, str_split_1(matched_peptides, ";"))
  }
  row <- tibble(ProteinId = protein_id)
  for (e in ENGINES) {
    current <- percolator_tibbles[[e]] %>% filter(ProteinId %in% ids_to_match)
    if (nrow(current) != 0) {
      col <- tibble(!!e := sum(current$num_peps))
    } else {
      col <- tibble(!!e := 0)
    }
    row <- bind_cols(row, col)
  }
  return(row)
}) %>%
  bind_rows()

# Determine which engines contributed the most to protein coverage
# Do so by multivariate linear regression -> it makes sense that
with_normal <- local({
  normal <- num_peptides_matched %>% select(all_of(c("ProteinId", names(normal_engine_tbs))))
  standard_search %>%
    select(., c(ProteinId, pcoverage_align)) %>%
    inner_join(., normal, by = join_by("ProteinId"))
})

reticulate::source_python(glue("{M$python_source}/trace_alignments.py"))
matched_peptides <- select(combined_results, c(ProteinId, MatchedPeptideIds))
protein_tb <- normal_engine_tbs %>%
  bind_rows() %>%
  separate_longer_delim(., "peptideIds", " ") %>%
  mutate(peptideIds = map_chr(peptideIds, clean_peptide))

measure_coverage_contribution <- function(protein_id) {
  row <- ENGINES
  return()
}



py$tracer <- AlignmentTracer(protein_tb, M$alignments$peptides, matched_peptides)
# Try to measure the proportion of alignments each engine takes up with its peptides
traced_intervals <- get_all_intervals(combined_results$ProteinId)
interval_coverage <- pmap(
  list(traced_intervals, combined_results$length, combined_results$ProteinId),
  \(x, y, z) interval_coverage(x, y, z)
) %>% bind_rows()
interval_groups <- pmap(
  list(traced_intervals, combined_results$ProteinId),
  \(x, y) get_interval_groups(x, y)
) %>%
  bind_rows() %>%
  as_tibble()

standard_engine_coverage <- interval_coverage %>%
  filter(ProteinId %in% standard_search$ProteinId) %>%
  select(-all_of(open_search_engines))

GRAPHS$engine_peptide_coverage <- standard_engine_coverage %>%
  pivot_longer(., cols = -ProteinId) %>%
  ggplot(aes(y = value, x = name, fill = name)) +
  geom_boxplot() +
  labs(title = "Density of per-engine peptide coverage") +
  ylab("Coverage (%)") +
  xlab("Engine name") +
  theme(legend.position = "none") +
  scale_fill_paletteer_d(PALETTE)

# Identify which engine, if any, is the best-performing
cov_list <- standard_engine_coverage %>%
  select(-ProteinId) %>%
  as.list()
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

#' Group engines that identify the same peptide groups using Jaccard distance
enginesXProtein <- peps <- num_peptides_matched %>%
  select(-all_of(open_search_engines)) %>%
  tb_transpose()

engine_dist <- vegan::vegdist(enginesXProtein, method = "jaccard")
dist_longer <- engine_dist |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = "engine") |>
  pivot_longer(cols = -engine)
GRAPHS$engine_sim_jaccard <- dist_longer |> ggplot(aes(x = engine, y = name, fill = value)) +
  geom_tile() +
  ylab("Engine") +
  xlab("Engine") +
  scale_fill_paletteer_c("ggthemes::Orange-Gold", name = "Jaccard distance")


#' Venn diagram for overlap
id_list <- lapply(normal_engine_tbs, \(x) x$ProteinId) %>%
  `names<-`(names(normal_engine_tbs))
venn <- ggVennDiagram(id_list, label = "none") +
  scale_fill_gradient(low = "grey90", high = "blue")
GRAPHS$engine_venn <- venn

save(c(GRAPHS, TABLES), glue("{M$outdir}/figures/engine_characteristics"))
