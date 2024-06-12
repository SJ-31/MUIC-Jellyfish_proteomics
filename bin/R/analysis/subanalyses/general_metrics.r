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
wanted <- c(
  "lineage", "Mods", "flashlfq_mean", "maxlfq_mean",
  "CAZy", "PFAMs", "EC", "BRITE",
  "interpro_accession", "PANTHER", "eggNOG_OGs", "UniProtKB_ID", "GO", "directlfq_mean"
)
GRAPHS$percent_found <- percent_found %>%
  dplyr::filter(!(first == 100 & sec == 100) & metric %in% wanted) %>%
  pivot_longer(cols = c("first", "sec")) %>%
  ggplot(aes(x = metric, y = value, fill = name)) +
  geom_bar(position = "dodge", stat = "identity") +
  ylab("% not missing") +
  scale_fill_discrete("Pass")

# Check if coverage and intensity differs significantly between protein groups
# for confirmation only (we expect them to differ)
tb <- data
grouping_metric <- "GO_category_MF"
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
  labs(y = "log intensity", x = grouping_metric)
GRAPHS$coverage_categories <- ggplotNumericDist(cov_list, "boxplot") +
  labs(y = "coverage (%)", x = grouping_metric)


with_category <- inner_join(tb, lfq) %>%
  select(ProteinId, log_intensity, !!grouping_metric) %>%
  filter(!is.na(log_intensity)) %>%
  arrange(log_intensity) %>%
  mutate(rank = seq_len(nrow(.)))
GRAPHS$category_ranks <- with_category %>%
  ggplot(aes(x = rank, y = log_intensity, color = !!grouping_metric)) +
  geom_point() +
  labs(x = "Rank", y = "Log intensity")

# Top ten most intense proteins
top_ten <- lfq %>%
  arrange(desc(log_intensity)) %>%
  slice(1:10)

# --------------------------------------------------------

## -# Annotation metrics
counts <- list()
counts$first <- getCounts(run$first)
counts$sec <- getCounts(run$sec)

wanted_cols <- c("GO_counts", "GO_max_sv", "num_peps", "pcoverage_align")
# Run Wilcox tests between on the metrics defined above, pairing up proteins that
# were identified in both runs
compare_tb <- inner_join(run$first, run$second, by = join_by(ProteinId), suffix = c(".first", ".sec"))

wilcox <- pairwiseFromTb(
  compare_tb, wanted_cols, c("less", "less", "less", "less"),
  \(x, y, ...) wilcox.test(x, y, paired = TRUE, na.rm = TRUE, ...)
)
sub <- substituteAll(
  c("num_peps", "pcoverage_align"),
  c("peptide number", "percent coverage"),
  \(x) gsub("_", " ", x)
)
TABLES$wilcox_pass_comparison <- wilcox %>%
  gt(rowname_col = "metric") %>%
  tab_header(
    title = md("**Wilcox test results**"),
    subtitle = "Comparison between first and second passes"
  ) %>%
  fmt(., columns = "metric", fns = \(x) map_chr(x, sub)) %>%
  text_case_match(
    "two.sided" ~ "two-sided",
    ".first" ~ "first",
    "_" ~ " pass ",
    .replace = "partial"
  ) %>%
  fmt_number(., columns = contains("p_value"), decimals = 5) %>%
  tab_stubhead(label = "Metric")

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
per_protein_wider <- local({
  p <- per_protein %>%
    pivot_longer(cols = c(first, sec), names_to = "pass")
  mean <- p %>% filter(type == "mean")
  stdev <- p %>% filter(type == "stdev")
  joined <- inner_join(mean, stdev, by = join_by(metric, pass)) %>%
    mutate(value = paste0(value.x, " (", value.y, ")")) %>%
    rename(percent_change = percent_change.x) %>%
    select(-matches("x|y"))
})
TABLES$per_protein_metrics <- per_protein_wider %>%
  pivot_wider(names_from = pass, values_from = value) %>%
  rename(
    second = sec,
    "Change from first to second (%)" = percent_change
  ) %>%
  gt() %>%
  tab_header(
    title = "Per protein metrics for each pass",
    subtitle = "Values reported as mean (stdev)"
  )

# Change in metrics from the first to the second pass
GRAPHS$per_protein_change <- per_protein %>%
  mutate(metric = paste0(metric, "_", type)) %>%
  pivot_longer(cols = c(first, sec)) %>%
  ggplot(aes(y = percent_change, x = metric, fill = metric)) +
  geom_bar(stat = "identity")

# Peptide characteristics
peptides <- getPeptideData(run$first$peptideIds)

# BUG Get rid of the largest peptides
peptides <- peptides |> filter(!length >= 500)
GRAPHS$peptide_lengths <- ggplot(peptides, aes(x = length)) +
  geom_histogram(binwidth = 10)

save(c(GRAPHS, TABLES), glue("{OUTDIR}/figures/general_metrics"))
