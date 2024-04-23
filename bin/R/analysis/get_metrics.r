library(gt)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)
args <- list(r_source = "./bin/R")
FDR <- 0.01
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))


EGGNOG_COLS <- c("EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")


PATH <- "./results/C_indra_A"
run <- runData("C_indra", TRUE, PATH)
OUTDIR <- glue("{PATH}/Analysis")
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR)
}
alignments <- alignmentData(PATH)
GRAPHS <- list()
TABLES <- list()


## -# Coverage metrics
cov_align <- compareFirstSecL(
  run, "pcoverage_align",
  TRUE, "ProteinId"
)
GRAPHS$run_coverage <- passDensityPlot(cov_align, 0.05) + labs(x = "percent coverage")


# Unique proteins to each run
run_uniques <- passUniques(run)
percent_found <- dplyr::bind_cols(
  notMissing(run$first),
  notMissing(run$sec)
) %>%
  `colnames<-`(c("first", "sec")) %>%
  tibble::rownames_to_column(., var = "metric") %>%
  as_tibble()
wanted <- c("lineage", "Mods", "flashlfq_mean", "maxlfq_mean", "CAZy", "PFAMs", "EC", "BRITE", "interpro_accession", "PANTHER", "eggNOG_OGs", "UniProtKB_ID", "GO", "directlfq_mean")
GRAPHS$percent_found <- percent_found %>%
  dplyr::filter(!(first == 100 & sec == 100) & metric %in% wanted) %>%
  pivot_longer(cols = c("first", "sec")) %>%
  ggplot(aes(x = metric, y = value, fill = name)) +
  geom_bar(position = "dodge", stat = "identity") +
  ylab("% not missing") +
  scale_fill_discrete("Pass")

# Check if coverage and intensity differs significantly between protein groups
# for confirmation only (we expect them to differ)
tb <- run$first
grouping_metric <- "category"
lfq <- dplyr::select(tb, all_of(grouping_metric), ProteinId) %>% inner_join(., mergeLfq(tb, "mean"))

apply_over <- tb[[grouping_metric]] %>%
  table() %>%
  discard(., \(x) x < 100) %>%
  names()
cov_list <- groupListFromTb(tb,
  v = apply_over, col_from = grouping_metric,
  target_col = "pcoverage_nmatch"
)
intensity_list <- groupListFromTb(lfq, apply_over, grouping_metric, "log_intensity")
GRAPHS$intensity_categories <- ggplotNumericDist(intensity_list, "boxplot") +
  labs(y = "log intensity", x = "category")
GRAPHS$coverage_categories <- ggplotNumericDist(cov_list, "boxplot") +
  labs(y = "coverage (%)", x = "category")


with_category <- inner_join(tb, lfq) %>%
  select(ProteinId, log_intensity, category) %>%
  filter(!is.na(log_intensity)) %>%
  arrange(log_intensity) %>%
  mutate(rank = seq_len(nrow(.)))
GRAPHS$category_ranks <- with_category %>%
  ggplot(aes(x = rank, y = log_intensity, color = category)) +
  geom_point() +
  labs(x = "Rank", y = "Log intensity")

# Top ten most intense proteins
top_ten <- lfq %>%
  arrange(desc(log_intensity)) %>%
  slice(1:10)


# --------------------------------------------------------

## PTM analysis
has_mods <- run$first %>% filter(!is.na(Mods))
UNIQUE_MODS <- run$first$Mods %>%
  discard(is.na) %>%
  map(\(x) str_split_1(x, "\\|")) %>%
  unlist() %>%
  map_chr(\(x) gsub(" [1-9]+", "", x)) %>%
  unique() %>%
  discard(\(x) grepl("0$", x))

ptms_first <- modMetrics(run$first)
ptms_sec <- modMetrics(run$sec)
ptm_percent_diff <- abs(ptms_sec$percentages - ptms_first$percentages)

# Test if modifications are associated using chi square

# Transform tb so that proteins are in rows and mods in cols
# A protein has TRUE in the mod col if it has that mod
tb <- run$first
chi <- tb %>%
  filter(category != "other") %>%
  select(ProteinId, category) %>%
  left_join(., df2Tb(ptms_first$count_df, "ProteinId")) %>%
  mutate(across(is.double, \(x) ifelse(is.na(x), FALSE, TRUE)))

no_mods <- chi %>%
  select(-c(ProteinId, category)) %>%
  apply(1, \(x) ifelse(any(x), FALSE, TRUE)) %>%
  unlist()
chi$none <- no_mods

# Create 2 x 2 contingency table with mods as columns and categories as rows
# But since mods are not mutually exclusive, will need to test each
# modification individually against the categories
# TODO: Need to format this nicely
chosen <- "Met_Oxidation"
tab <- table(chi$category, chi[[chosen]])
ptm_tests <- list()
for (mod in UNIQUE_MODS) {
  ptm_tests[[mod]]$test <- chisq.test(chi$category, chi[[mod]])
  ptm_tests[[mod]]$table <- table(chi$category, chi[[mod]])
  print(ptm_tests[[mod]])
}

# --------------------------------------------------------


## -# Annotation metrics
counts <- list()
counts$first <- getCounts(run$first)
counts$sec <- getCounts(run$sec)

avgStdevs <- function(tb, cols, stat) {
  result <- vector()
  for (col in cols) {
    result <- c(result, stat(tb[[col]]))
  }
  return(result)
}


wanted_cols <- c("GO_counts", "GO_max_sv", "num_peps", "pcoverage_align")
# Run Wilcox tests between on the metrics defined above, pairing up proteins that
# were identified in both runs
WILCOX <- lapply(wanted_cols, \(x) {
  compare_tb <- compareFirstSecW(run, x, "ProteinId", TRUE)
  test_tb <- wilcoxWrapper(compare_tb, TRUE, x)
  return(test_tb)
}) %>%
  bind_rows() %>%
  mutate(alternative = map_chr(alternative, \(x)
  ifelse(x == "two.sided", x, paste0("first_", x))))


# Results per protein
# Evaluate significance of each
per_protein <- tibble(
  metric = rep(wanted_cols, 2),
  type = c(rep("mean", length(wanted_cols)), rep("stdev", length(wanted_cols))),
  first = c(
    avgStdevs(run$first, wanted_cols, \(x) mean(x, na.rm = TRUE)),
    avgStdevs(run$first, wanted_cols, \(x) sd(x, na.rm = TRUE))
  ),
  sec = c(
    avgStdevs(run$sec, wanted_cols, \(x) mean(x, na.rm = TRUE)),
    avgStdevs(run$sec, wanted_cols, \(x) sd(x, na.rm = TRUE))
  )
) %>%
  mutate(percent_change = (sec - first) / first) %>%
  mutate(across(is.double, \(x) round(x, 3)))
per_protein %>%
  gt() %>%
  tab_header("Comparison between first and second passes") %>%
  gtsave(glue("{OUTDIR}/per_protein_metrics.tex"))
write_tsv(per_protein, file = glue("{OUTDIR}/per_protein_metrics.tsv"))

# Change in metrics from the first to the second pass
GRAPHS$per_protein_change <- per_protein %>%
  mutate(metric = paste0(metric, "_", type)) %>%
  pivot_longer(cols = c(first, sec)) %>%
  ggplot(aes(y = percent_change, x = metric, fill = metric)) +
  geom_bar(stat = "identity")


# Compare with previous results
p_rename <- c(NCBI_ID = "Accession Number", pcoverage_nmatch.prev = "Sequence coverage [%]")
p_all <- read_tsv("./data/reference/previous_all.tsv") %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))
p_toxins <- read_tsv("./data/reference/previous_toxins.tsv") %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))

compare_all <- inner_join(p_all, run$first, by = join_by(NCBI_ID)) %>%
  left_join(., run$sec, by = join_by(NCBI_ID), suffix = JOIN_SUFFIX) %>%
  select(c(NCBI_ID, header.first, category.first, pcoverage_nmatch.prev, pcoverage_nmatch.first, pcoverage_nmatch.sec)) %>%
  mutate(
    pcoverage_nmatch.first = pcoverage_nmatch.first * 100,
    pcoverage_nmatch.sec = pcoverage_nmatch.sec * 100
  )

cov_longer <- compare_all %>%
  select(contains("coverage"), NCBI_ID) %>%
  rename_with(., \(x) purrr::map_chr(x, \(y) {
    if (grepl("NCBI_ID", y)) {
      return(y)
    }
    y <- gsub("pcoverage_nmatch", "", y)
    if (y == ".prev") {
      return("previous")
    } else if (y == ".first") {
      return("first")
    }
    return("second")
  })) %>%
  mutate(
    first_diff = first - previous,
    second_diff = second - previous
  ) %>%
  select(-c(first, second, previous)) %>%
  pivot_longer(cols = !NCBI_ID) %>%
  mutate(value = round(value, 2))

GRAPHS$shared_cov_bp <- cov_longer %>% ggplot(aes(y = value, fill = name)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  xlab("Seqence coverage") +
  ggtitle("Shared protein sequence coverage")

# For space purposes, drop entries where ALL had less than 30% sequence coverage, as well as any entries both second and first passes cannot be compared
threshold <- 0
cov_longer <- cov_longer %>%
  group_by(NCBI_ID) %>%
  nest()
mask <- cov_longer %>% apply(1, \(x) {
  data <- x$data
  if (all(data$value < threshold) | any(is.na(data$value))) {
    return(FALSE)
  }
  return(TRUE)
})
cov_longer <- cov_longer[mask, ] %>% unnest()
cov_longer$name <- factor(cov_longer$name, levels = c("previous", "first", "second"))

GRAPHS$protein_wise_coverage <- cov_longer %>%
  ggplot(aes(x = NCBI_ID, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Sequence coverage change from previous (%)") +
  xlab(element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )

p_test <- wilcox.test(compare_all$pcoverage_nmatch.prev, compare_all$pcoverage_nmatch.first, paired = TRUE, alternative = "less")
WILCOX <- bind_rows(WILCOX, tibble(
  alternative = "prev_less_than_first",
  metric = "coverage",
  statistic = p_test$statistic,
  p_value = p_test$p.value,
  reject_null = ifelse(p_test$p.value < 0.05, "Y", "N")
))
write_tsv(WILCOX, file = glue("{OUTDIR}/wilcox_tests.tsv"))


# Filter out new proteins
compare_toxin <- compare_all %>% filter(NCBI_ID %in% p_toxins$NCBI_ID)
new_proteins <- run$first %>% filter(!NCBI_ID %in% compare_all$NCBI_ID)
new_toxins <- new_proteins %>% filter(category == "venom_component")

# --------------------------------------------------------
# Lineage metrics


# --------------------------------------------------------

num_peps <- compareFirstSecW(run, "num_peps", "ProteinId", TRUE)

## -#  Engine analysis
##
tb <- run$first
num_ids <- tb %>%
  filter(ProteinGroupId != "U") %>%
  select(c(ProteinGroupId, pcoverage_nmatch, ProteinId, num_peps, num_unique_peps)) %>%
  mutate(engine_count = purrr::map_dbl(ProteinGroupId, \(x) {
    x <- splitGroupStr(x, TRUE, TRUE) %>%
      discard(\(x) x == "U")
    return(length(x))
  }))

# Figure out which engines had the biggest contributions
engine_counts <- num_ids$ProteinGroupId %>%
  lapply(., \(x) splitGroupStr(x, TRUE)) %>%
  unlist() %>%
  discard(\(x) x == "U") %>%
  table()

# Is there an association between the number of matched peptides and the
# identity of the engines matching them? That is, are some engines more
# "isolated" than others? i.e. engine B tends to identify proteins that other
# engine don't. But also remember that in your setup each protein needs to
# be identified by at least two standard engines (only one in open search)
tb <- run$first
engine_hits <- tb %>%
  select(ProteinId, ProteinGroupId) %>%
  apply(1, \(x) {
    groups <- splitGroupStr(
      x["ProteinGroupId"],
      TRUE, TRUE
    )
    row <- tibble(ProteinId = x["ProteinId"])
    for (g in groups) {
      row <- add_column(row, !!g := TRUE)
    }
    return(row)
  }) %>%
  bind_rows() %>%
  replaceNaAll(., FALSE)
# Collapse the contingency table?
engine_hits <- engine_hits %>%
  inner_join(., dplyr::select(tb, category, ProteinId)) %>%
  select(-ProteinId) %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  filter(n() > 30) %>%
  ungroup()

categories <- unique(engine_hits$category)

# An example of what the contingency table for the tests would look like,
# and what the expected values would be under the chi-square test of independence
current <- "comet"
tl <- table(engine_hits$category, engine_hits[[current]])
ct <- "membrane"
show_contigency <- table(engine_hits$category == ct, engine_hits[[current]])
expected <- formatEngineContingency(show_contigency, ct, TRUE)
expected %>%
  gt() %>%
  tab_header(
    title = glue("Engine: {current}, Category: {ct}"),
    subtitle = "Expected values under chi-square in parentheses"
  )

engines <- colnames(engine_hits)[!colnames(engine_hits) %in% c("U", "category")]
engine_chi <- tibble()
engine_table_list <- list()
for (engine in engines) {
  engine_table_list[[engine]] <- list()
  for (cat in categories) {
    t <- table(engine_hits$category == cat, engine_hits[[engine]])
    engine_table_list[[engine]][[cat]] <- t
    print(formatEngineContingency(t, TRUE))
    tst <- chisq.test(t)
    row <- tibble(engine = engine, category = cat, p_value = tst$p.value)
    engine_chi <- bind_rows(engine_chi, row)
  }
}
engine_chi <- engine_chi %>% mutate(
  p_adjust = p.adjust(p_value, "holm"),
  is_significant = ifelse(p_value < 0.05, TRUE, FALSE)
)

# Compute effect size for significant hits, using the odds ratio
# This represents how much more likely the engine is to identify
# peptides from a protein of the given category
odds_ratios <- engine_chi %>%
  apply(1, \(x) {
    if (!as.logical(x["is_significant"])) {
      return(tibble(OR = NA, OR_upper_ci = NA, OR_lower_ci = NA))
    }
    e <- x["engine"]
    c <- x["category"]
    table <- engine_table_list[[e]][[c]]
    odds_ratio <- oddsRatio(table)
    odds_ratio_ci_u <- oddsRatio(table, TRUE, "upper")
    odds_ratio_ci_l <- oddsRatio(table, TRUE, "lower")
    row <- tibble(
      OR = odds_ratio,
      OR_upper_ci = odds_ratio_ci_u,
      OR_lower_ci = odds_ratio_ci_l
    )
    return(row)
  }) %>%
  bind_rows()
engine_chi <- bind_cols(engine_chi, odds_ratios)
engine_chi %>%
  gt() %>%
  tab_header(
    title = "Association between engine hits and protein GO category",
    subtitle = "Tested using chi-square"
  ) %>%
  gtsave(., glue("{OUTDIR}/engine_association_results.html"))


# Record all results
engine_table_tb <- tibble(category = categories)
remove_cat <- FALSE
for (e in names(engine_table_list)) {
  e_list <- engine_table_list[[e]]
  temp_lst <- list()
  remove_col_label <- FALSE
  for (ct in names(e_list)) {
    tb <- formatEngineContingency(e_list[[ct]], w_expected = TRUE)
    if (remove_cat) {
      tb <- dplyr::select(tb, -category)
    }
    tb <- gt(tb)
    if (!remove_col_label) {
      remove_col_label <- TRUE
    } else {
      tb <- tb %>% tab_options(column_labels.hidden = TRUE)
    }
    temp_lst <- c(temp_lst, as_raw_html(tb))
  }
  temp_tb <- tibble(!!e := temp_lst)
  engine_table_tb <- bind_cols(engine_table_tb, temp_tb)
  remove_cat <- TRUE
}
engine_table_tb %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  tab_header(
    title = md("**Contingency tables used in chi-square analysis of engine bias**"),
    subtitle = "For each table, misses are on the left, hits on the right.
Expected values are in parentheses"
  ) %>%
  gtsave(., glue("{OUTDIR}/engine_category_contingency.html"))


# Correlation between no. identifications by engines and coverage
engine_cor <- cor.test(num_ids$engine_count, num_ids$pcoverage_nmatch) # A weak positive correlation, but statistically significant
# Correlation between peptide number and coverage
n_peps_cor <- cor.test(num_ids$num_unique_peps, num_ids$pcoverage_nmatch) # A weak positive correlation, but statistically significant


# --------------------------------------------------------

# Grouping proteins by their participation in KEGG pathways
current <- run$first
pwy <- groupPathways(current)

# --------------------------------------------------------

# --------------------------------------------------------
# Investigating trends in missing quantification
current <- run$first
cq <- c("directlfq", "maxlfq", "flashlfq")
missing_quant_tests <- list()
for (q in cq) {
  noq <- current %>% filter(is.na(!!as.symbol(glue("{cq}_mean"))))
  hasq <- current %>% filter(!ProteinId %in% noq$ProteinId)
  missing_quant_tests[[q]] <- wilcox.test(noq$num_peps, hasq$num_peps,
    alternative = "l"
  )
}
capture.output(missing_quant_tests, file = glue("{OUTDIR}/missing_quantification_tests.txt"))
# --------------------------------------------------------
