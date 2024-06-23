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
  runData(runs[[x]],
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

getCommonColumn <- function(tb_list, target_col, predicate) {
  lapply(tb_list, \(x) {
    x <- predicate(x)
    x[[target_col]]
  })
}


common_ids <- accumulate(id_list, intersect)[[length(id_list)]]

coverage_comparison <- getCommonColumn(
  data, "pcoverage_align",
  \(x) dplyr::filter(x, x$ProteinId %in% common_ids & !is.na(pcoverage_align))
)

num_peps <- getCommonColumn(
  data, "num_peps",
  \(x) dplyr::filter(x, x$ProteinId %in% common_ids & !is.na(num_peps))
)

num_unique_peps <- getCommonColumn(
  data, "num_unique_peps",
  \(x) dplyr::filter(x, x$ProteinId %in% common_ids & !is.na(num_unique_peps))
)

testAndVisualize <- function(data_list, plot_type, pairwiseTest) {
  results <- list(graph = NULL, pairwise = NULL, kruskal = NULL)
  results$kruskal <- kruskal.test(data_list)
  if (results$kruskal$p.value < 0.05) {
    results$pairwise <- testAllPairs(data_list, pairwiseTest) %>% mutate(
      p_adjust = p.adjust(p_value),
      is_significant = ifelse(p_adjust < 0.05, "Y", "N")
    )
  }
  results$graph <- ggplotNumericDist(data_list, plot_type)
  return(results)
}

cc <- testAndVisualize(coverage_comparison, "boxplot", wilcox.test)

coverage_comparison_g <- cc$graph + ylab("coverage (%)") +
  scale_color_paletteer_d(PALETTE) +
  guides(color = guide_legend("Run parameter")) +
  theme(
    axis.text.x = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 22), legend.title = element_text(size = 20, face = "bold")
  )

legend <- getLegend(coverage_comparison_g)

coverage_comparison_g <- coverage_comparison_g + guides(color = "none")

peptide_number <- testAndVisualize(num_peps, "boxplot", wilcox.test)
peptide_number_g <- peptide_number$graph +
  ylab("Matched peptide number per protein") +
  scale_color_paletteer_d(PALETTE) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  guides(color = "none")

unique_peptide_number <- testAndVisualize(
  num_unique_peps, "boxplot", wilcox.test
)
unique_peptide_number_g <- unique_peptide_number$graph +
  ylab("Matched unique peptide number per protein") + guides(color = "none") +
  scale_color_paletteer_d(PALETTE) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + guides(color = "none")


TABLES$all_params_tests <- bind_rows(
  cc$pairwise |> mutate(metric = "coverage"),
  peptide_number$pairwise |> mutate(metric = "peptide number"),
  unique_peptide_number$pairwise |> mutate(metric = "unique peptide number")
) |>
  gt()


GRAPHS$coverage_metrics <- grid.arrange(legend, coverage_comparison_g,
  peptide_number_g, unique_peptide_number_g,
  nrow = 2, ncol = 2
)

#' ---------------------------------------------------
# Analyze GO and KEGG pathway composition
GO_list <- lapply(data, \(x) {
  list(GO_IDs = goVector(x, go_column = "GO_IDs"), GO_slims = goVector(x, go_column = "GO_slims"))
})
unique_terms <- GO_list %>%
  unlist() %>%
  unique()

go_slims <- local({
  tbs <- lmap(GO_list, \(x) {
    table2Tb(table(x[[1]]$GO_slims), "GO_IDs") %>%
      rename(!!as.symbol(names(x)) := n)
  })
  purrr::reduce(tbs, \(acc, x) left_join(acc, x), .init = tibble(GO_IDs = unique_terms))
}) %>%
  filter(!if_all(where(is.numeric), is.na)) %>%
  mutate(across(is.numeric, \(x) replace(x, is.na(x), 0)))

paramXSlim <- tbTranspose(go_slims)
dist <- vegan::vegdist(paramXSlim, method = "bray")
pcoa <- vegan::wcmdscale(dist, eig = TRUE)
GRAPHS$param_pcoa <- pcoa$points |>
  as.data.frame() |>
  rownames_to_column(var = "Parameter") |>
  ggplot(aes(x = Dim1, y = Dim2, color = Parameter)) +
  geom_point() +
  scale_color_paletteer_d(PALETTE) +
  xlab("V1") +
  ylab("V2")

go_terms <- local({
  tbs <- lmap(GO_list, \(x) {
    table2Tb(table(x[[1]]$GO_IDs), "GO_IDs") %>%
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

compare_cols <- c("pcoverage_nmatch", "pcoverage_align", "num_unique_peps", "num_peps")
merged <- merged %>% mutate(across(contains(compare_cols), as.numeric))
tests <- pairwiseFromTb(
  merged,
  c("pcoverage_nmatch", "pcoverage_align", "num_unique_peps", "num_peps"),
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
  as.list()
names(peps_list) <- replaceSuffix(
  names(peps_list),
  c(
    ".def" = " default", ".nd" = " no denovo",
    "num_peps" = "", "num_unique_peps" = "unique"
  )
)


GRAPHS$default_no_denovo_coverage <- ggplotNumericDist(coverage_list, method = "freq_poly", linewidth = 0.9) +
  xlab("Alignment coverage (%)") +
  ylab("Frequency") +
  labs(
    title = "Coverage comparison between runs with and without added de novo peptides",
    color = "Parameter"
  ) + scale_color_paletteer_d(PALETTE2)


GRAPHS$default_no_denovo_peptide_count <- ggplotNumericDist(peps_list, method = "boxplot") +
  xlab("Run parameter") +
  ylab("Count") +
  labs(
    title = "Number of identified peptides between runs with and without added de novo peptides",
    color = "Parameter"
  ) + scale_color_paletteer_d(PALETTE2)

TABLES$denovo_test <- tests |>
  mutate(alternative = map_chr(alternative, \(x) str_replace(x, "\\.def", "default"))) |>
  gt()

save(c(GRAPHS, TABLES), glue("{M$outdir}/figures/run_parameters"))
