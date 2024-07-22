if (!exists("SOURCED")) {
  source(paste0(getwd(), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
PALETTE <- "ggthemes::hc_default"
PALETTE2 <- "ggthemes::Classic_Color_Blind"

library("tidyverse")
library("gridExtra")
library("ggVennDiagram")
library("ggplot2")
library("ggridges")
library("venn")
library("Peptides")
library("glue")

CHOSEN_PASS <- "first"
args <- list(
  r_source = glue("{M$wd}/bin/R"),
  python_source = glue("{M$wd}/bin")
)
source(glue("{M$r_source}/helpers.r"))
source(glue("{M$r_source}/GO_helpers.r"))
source(glue("{M$r_source}/cluster_helpers.r"))
source(glue("{M$r_source}/DR_helpers.r"))
source(glue("{M$r_source}/analysis/metric_functions.r"))


runs <- list("C_indra.calibrated", "C_indra.msconvert", "C_indra", "ND_C_indra")
names(runs) <- c("msgf_calibrated", "peak_picked", "default", "no_denovo")
data <- lapply(names(runs), \(x) {
  get_run(runs[[x]],
    glue("{M$wd}/results/{runs[[x]]}"),
    which = CHOSEN_PASS
  )
}) %>% `names<-`(names(runs))

# Decide which run parameters are best

GRAPHS <- list()
TABLES <- list()

# Intersection of protein_ids
id_list <- lapply(data, \(x) x$ProteinId)
venn <- ggVennDiagram(id_list, label_alpha = 0) +
  scale_fill_gradient(low = "#eff1f5", high = "#40a02b")
GRAPHS$run_param_venn <- venn

get_common_col <- function(tb_list, target_col, predicate) {
  lapply(tb_list, \(x) {
    x <- predicate(x)
    x[[target_col]]
  })
}


common_ids <- accumulate(id_list, intersect)[[length(id_list)]]

coverage_comparison <- get_common_col(
  data, "pcoverage_align",
  \(x) dplyr::filter(x, x$ProteinId %in% common_ids & !is.na(pcoverage_align))
)

num_peps <- get_common_col(
  data, "num_peps",
  \(x) dplyr::filter(x, x$ProteinId %in% common_ids & !is.na(num_peps))
)

num_unique_peps <- get_common_col(
  data, "num_unique_peps",
  \(x) dplyr::filter(x, x$ProteinId %in% common_ids & !is.na(num_unique_peps))
)

test_and_visualize <- function(data_list, plot_type, transform_fn = NULL) {
  results <- list(graph = NULL, pairwise = NULL, kruskal = NULL, pairwise_conclusion = NULL)
  results$kruskal <- kruskal.test(data_list)
  if (results$kruskal$p.value < 0.05) {
    results$pairwise <- test_all_pairs(
      data_list,
      \(x, y) wilcox.test(x, y, alternative = "greater"), "greater"
    ) |>
      bind_rows(test_all_pairs(data_list, wilcox.test, two_sided = TRUE)) |>
      mutate(p_adjust = p.adjust(p_value), significant = map_dbl(p_adjust, \(x) x < 0.05))
    results$pairwise_conclusion <- conclude_one_sided(results$pairwise)
  }
  if (!is.null(transform_fn)) {
    data_list <- transform_fn(data_list)
  }
  results$graph <- gg_numeric_dist(data_list, plot_type)
  return(results)
}


cc <- test_and_visualize(coverage_comparison, "boxplot", transform_fn = \(x) lapply(x, \(y) log(y)))

TABLES$coverage_comparison_conclusion <- gt(cc$pairwise_conclusion)
TABLES$coverage_comparison <- cc$pairwise

coverage_comparison_g <- cc$graph + ylab("log Percent coverage") +
  scale_color_paletteer_d(PALETTE) +
  guides(color = guide_legend("Run parameter")) +
  theme(
    axis.text.x = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 22), legend.title = element_text(size = 20, face = "bold")
  )

legend <- get_legend(coverage_comparison_g)

coverage_comparison_g <- coverage_comparison_g + guides(color = "none")

peptide_number <- test_and_visualize(num_peps, "boxplot",
  transform_fn = \(x) lapply(x, \(y) log(y))
)
peptide_number_g <- peptide_number$graph +
  ylab("log Matched peptide number per protein") +
  scale_color_paletteer_d(PALETTE) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  guides(color = "none")

TABLES$peptide_number_conclusion <- gt(peptide_number$pairwise_conclusion)
TABLES$peptide_number <- peptide_number$pairwise

unique_peptide_number <- test_and_visualize(
  num_unique_peps, "boxplot",
  transform_fn = \(x) lapply(x, \(y) log(y))
)
unique_peptide_number_g <- unique_peptide_number$graph +
  ylab("log Matched unique peptide number per protein") + guides(color = "none") +
  scale_color_paletteer_d(PALETTE) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + guides(color = "none")

TABLES$unique_peptide_number_conclusion <- gt(unique_peptide_number$pairwise_conclusion)
TABLES$unique_peptide_number <- unique_peptide_number$pairwise

GRAPHS$coverage_metrics <- grid.arrange(legend, coverage_comparison_g,
  peptide_number_g, unique_peptide_number_g,
  nrow = 2, ncol = 2
)

#' ---------------------------------------------------
# Analyze GO and KEGG pathway composition
GO_list <- lapply(data, \(x) {
  list(GO_IDs = get_go_vec(x, go_column = "GO_IDs"), GO_slims = get_go_vec(x, go_column = "GO_slims"))
})
unique_terms <- GO_list %>%
  unlist() %>%
  unique()

go_slims <- local({
  tbs <- lmap(GO_list, \(x) {
    table2tb(table(x[[1]]$GO_slims), "GO_IDs") %>%
      rename(!!as.symbol(names(x)) := n)
  })
  purrr::reduce(tbs, \(acc, x) left_join(acc, x), .init = tibble(GO_IDs = unique_terms))
}) %>%
  filter(!if_all(where(is.numeric), is.na)) %>%
  mutate(across(is.numeric, \(x) replace(x, is.na(x), 0)))

paramXSlim <- tb_transpose(go_slims)
dist <- vegan::vegdist(paramXSlim, method = "bray")
dist_longer <- dist |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = "parameter") |>
  pivot_longer(cols = -parameter)
GRAPHS$param_heatmap <- dist_longer |> ggplot(aes(x = parameter, y = name, fill = value)) +
  geom_tile() +
  ylab("Parameter") +
  xlab("Parameter") +
  scale_fill_paletteer_c("ggthemes::Purple", name = "Bray-Curtis distance")

go_terms <- local({
  tbs <- lmap(GO_list, \(x) {
    table2tb(table(x[[1]]$GO_IDs), "GO_IDs") %>%
      rename(!!as.symbol(names(x)) := n)
  })
  purrr::reduce(tbs, \(acc, x) left_join(acc, x), .init = tibble(GO_IDs = unique_terms))
}) %>%
  mutate(across(is.numeric, \(x) replace(x, is.na(x), 0)))

#' ----------
#' Evaluate contribution of denovo peptides
merged <- inner_join(data$default, data$no_denovo,
  by = join_by(ProteinId),
  suffix = c(".def", ".nd")
)
coverage <- list(default = merged$pcoverage_align.def, no_denovo = merged$pcoverage_align.nd)

compare_cols <- c("pcoverage_align", "pcoverage_align", "num_unique_peps", "num_peps")
merged <- merged %>% mutate(across(contains(compare_cols), as.numeric))
tests <- pairwise_tests_tb(
  merged,
  c("pcoverage_align", "pcoverage_align", "num_unique_peps", "num_peps"),
  rep("greater", 4),
  \(x, y, ...) wilcox.test(x, y, paired = TRUE, ...),
  suffixes = c(".def", ".nd")
)

replaceSuffix <- function(name_list, replacement_map) {
  map_chr(name_list, \(x) str_replace_all(x, replacement_map))
}

coverage_list <- merged %>%
  select(contains("coverage_align")) %>%
  as.list()
names(coverage_list) <- replaceSuffix(
  names(coverage_list),
  c(".def" = " default", ".nd" = " no denovo", "pcoverage_align" = "")
)

peps_list <- merged %>%
  select(contains("peps")) %>%
  as.list() |>
  lapply(\(x) log(x))
names(peps_list) <- replaceSuffix(
  names(peps_list),
  c(
    ".def" = " default", ".nd" = " no denovo",
    "num_peps" = "", "num_unique_peps" = "unique"
  )
)


GRAPHS$default_no_denovo_coverage <- gg_numeric_dist(coverage_list, method = "freq_poly", linewidth = 0.9) +
  xlab("Alignment coverage (%)") +
  ylab("Frequency") +
  labs(
    title = "Coverage comparison between runs with and without added de novo peptides",
    color = "Parameter"
  ) + scale_color_paletteer_d(PALETTE2)


GRAPHS$default_no_denovo_coverage_2 <- ggplot(
  merged,
  aes(x = pcoverage_align.def, y = pcoverage_align.nd, color = length.nd)
) +
  geom_point() +
  paletteer::scale_colour_paletteer_c("grDevices::Cold") +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
    linetype = 2,
    colour = "black"
  ) +
  ylab("Coverage without de novo peptides") +
  xlab("Coverage with de novo peptides") +
  guides(color = guide_legend("Protein length"))


GRAPHS$default_no_denovo_peptide_count <- gg_numeric_dist(peps_list, method = "boxplot") +
  xlab("Run parameter") +
  ylab("log Count") +
  labs(
    title = "Number of identified peptides between runs with and without added de novo peptides",
    color = "Parameter"
  ) + scale_color_paletteer_d(PALETTE2)

TABLES$denovo_test <- tests |>
  mutate(alternative = map_chr(alternative, \(x) str_replace(x, "\\.def", "default"))) |>
  gt()

save(c(GRAPHS, TABLES), glue("{M$outdir}/Figures/run_parameters"))
