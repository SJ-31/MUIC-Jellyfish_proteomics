library("tidyverse")
library("ggVennDiagram")
library("ggplot2")
library("ggridges")
library("venn")
library("Peptides")
library("glue")
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}

CHOSEN_PASS <- "first"
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  embd_type = "protein",
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo")
)
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))

FDR <- 0.05
PATH <- glue("{wd}/results/")

runs <- list("C_indra.calibrated", "C_indra.msconvert", "C_indra", "ND_C_indra")
names(runs) <- c("msgf_calibrated", "peak_picked", "default", "no_denovo")
data <- lapply(names(runs), \(x) {
  runData(runs[[x]],
    glue("{PATH}/{runs[[x]]}"),
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
GRAPHS$coverage_comparison <- cc$graph + ylab("coverage (%)") + xlab("Run modification")

peptide_number <- testAndVisualize(num_peps, "boxplot", wilcox.test)
GRAPHS$peptide_number <- peptide_number$graph +
  ylab("Matched peptide number per protein")

unique_peptide_number <- testAndVisualize(
  num_unique_peps, "boxplot", wilcox.test
)
GRAPHS$unique_peptide_number <- unique_peptide_number$graph +
  ylab("Matched unique peptide number per protein")


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
# pcoa$points %>% ggplot(aes(x = ))

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
  select(contains("coverage")) %>%
  as.list()
names(coverage_list) <- replaceSuffix(
  names(coverage_list),
  c(".def" = " default", ".nd" = " no denovo")
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
  )
GRAPHS$default_no_denovo_coverage

GRAPHS$default_no_denovo_peptide_count <- ggplotNumericDist(peps_list, method = "boxplot") +
  xlab("Run parameter") +
  ylab("Count") +
  labs(
    title = "Number of identified peptides between runs with and without added de novo peptides",
    color = "Parameter"
  )
