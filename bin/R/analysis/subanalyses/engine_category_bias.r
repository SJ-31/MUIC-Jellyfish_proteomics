if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

GRAPHS <- list()
TABLES <- list()

#'   Engine category bias
#'
#' 1) Is there an association between the number of matched peptides and the
#' identity of the engines matching them? That is, are some engines more
#' "isolated" than others? i.e. engine B tends to identify proteins that other
#' engine don't. But also remember that in your setup each protein needs to
#' be identified by at least two standard engines (only one in open search)


num_peps <- compareFirstSecW(M$run, "num_peps", "ProteinId", TRUE)
tb <- M$data
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

compare_col <- "GO_category_MF"
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
hits <- engine_hits

# Collapse the contingency table
engine_hits <- engine_hits %>%
  inner_join(., dplyr::select(tb, !!compare_col, ProteinId)) %>%
  select(-ProteinId) %>%
  filter(!is.na(!!compare_col)) %>%
  group_by(!!as.symbol(compare_col)) %>%
  filter(n() > 30) %>%
  ungroup()

categories <- unique(engine_hits[[compare_col]])

# An example of what the contingency table for the tests would look like,
# and what the expected values would be under the chi-square test of independence
current <- "comet"
tl <- table(engine_hits[[compare_col]], engine_hits[[current]])
ct <- "protein binding"
show_contigency <- table(engine_hits[[compare_col]] == ct, engine_hits[[current]])
expected <- formatEngineContingency(show_contigency, ct, TRUE)
expected %>%
  gt() %>%
  tab_header(
    title = glue("Engine: {current}, Category: {ct}"),
    subtitle = "Expected values under chi-square in parentheses"
  )
TABLES$sample_expected <- expected

engines <- colnames(engine_hits)[!colnames(engine_hits) %in% c("U", compare_col)]
engine_chi <- tibble()
engine_table_list <- list()
for (engine in engines) {
  engine_table_list[[engine]] <- list()
  for (cat in categories) {
    t <- table(engine_hits[[compare_col]] == cat, engine_hits[[engine]])
    engine_table_list[[engine]][[cat]] <- t
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
# Confidence interval is 95%
odds_ratios <- engine_chi %>%
  apply(1, \(x) {
    row <- tibble(OR = NA, OR_upper_ci = NA, OR_lower_ci = NA)
    if (!as.logical(x["is_significant"])) {
      return(row)
    }
    e <- x["engine"]
    c <- x["category"]
    table <- engine_table_list[[e]][[c]]
    try({
      odds_ratio <- oddsRatio(table)
      odds_ratio_ci_u <- oddsRatio(table, TRUE, "upper")
      odds_ratio_ci_l <- oddsRatio(table, TRUE, "lower")
      row <- tibble(
        OR = odds_ratio,
        OR_upper_ci = odds_ratio_ci_u,
        OR_lower_ci = odds_ratio_ci_l
      )
    })
    return(row)
  }) %>%
  bind_rows()
engine_chi_f <- bind_cols(engine_chi, odds_ratios)
TABLES$engine_chi <- engine_chi_f %>%
  mutate(across(is.double, \(x) round(x, 4))) %>%
  mutate(
    OR = paste0(OR, " [", OR_lower_ci, ", ", OR_upper_ci, "]"),
    OR = map_chr(OR, \(x) ifelse(grepl("NA|NaN", x), "NA", x))
  ) %>%
  select(-matches("OR_")) %>%
  rename("Odds ratio, 95% CI [lower, upper]" = OR) %>%
  gt() %>%
  tab_header(
    title = "Association between engine hits and protein GO category",
    subtitle = "Tested using chi-square"
  )
# Odds ratio is interpreted as the given engine is OR times as likely to miss
# the protein when it is NOT of the stated category than when it is
# So higher OR means that the engine is biased towards hitting proteins in that category


# Record all results
engine_table_tb <- tibble(category = categories)
remove_cat <- FALSE
for (e in names(engine_table_list)) {
  e_list <- engine_table_list[[e]]
  temp_lst <- list()
  remove_col_label <- FALSE
  for (ct in names(e_list)) {
    tib <- formatEngineContingency(e_list[[ct]], w_expected = TRUE)
    if (remove_cat) {
      tib <- dplyr::select(tib, -category)
    }
    tib <- gt(tib)
    if (!remove_col_label) {
      remove_col_label <- TRUE
    } else {
      tib <- tib %>% tab_options(column_labels.hidden = TRUE)
    }
    temp_lst <- c(temp_lst, as_raw_html(tib))
  }
  temp_tb <- tibble(!!e := temp_lst)
  engine_table_tb <- bind_cols(engine_table_tb, temp_tb)
  remove_cat <- TRUE
}
TABLES$engine_category_contingency <- engine_table_tb %>%
  gt() %>%
  fmt_markdown(columns = everything()) %>%
  tab_header(
    title = md("**Contingency tables used in chi-square analysis of engine bias**"),
    subtitle = "For each table, misses are on the left, hits on the right.
Expected values are in parentheses"
  )


# Correlation between no. identifications by engines and coverage
engine_cor <- cor.test(num_ids$engine_count, num_ids$pcoverage_nmatch)
engine_cor$data.name <- "Correlation between number of identifications by different engines and coverage"
# A weak positive correlation, but statistically significant
# Correlation between peptide number and coverage
n_peps_cor <- cor.test(num_ids$num_unique_peps, num_ids$pcoverage_nmatch)
engine_cor$data.name <- "Correlation between number of identified peptides and coverage"
TABLES$correlation <- gt(bind_rows(htest2Tb(engine_cor), htest2Tb(n_peps_cor)))
# A weak positive correlation, but statistically significant



intensity <- mergeLfq(tb, "mean") %>%
  filter(!is.na(log_intensity)) |>
  select(ProteinId, log_intensity)
intensity$intensity_class <- map_chr(intensity$log_intensity, \(x) {
  if (x <= quantile(intensity$log_intensity, 0.25)) {
    "first"
  } else if (quantile(intensity$log_intensity, 0.75) <= x) {
    "third"
  } else {
    "second"
  }
})
hits <- hits |> inner_join(intensity, by = join_by(ProteinId))

intense <- chisqNME(
  tb = hits, var_a_levels = engines,
  var_b_col = "intensity_class", var_a = "engine", var_b = "intensity_class"
)

TABLES$engine_intensity_chi <- intense$gt$chi
TABLES$engine_intensity_contingency <- intense$gt$contingency

save(c(TABLES, GRAPHS), glue("{M$outdir}/figures/engine_category_bias"))
