if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
ta <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = ta)

GRAPHS <- list()
TABLES <- list()

#'   Engine category bias
#'
#' 1) Is there an association between the number of matched peptides and the
#' identity of the engines matching them? That is, are some engines more
#' "isolated" than others? i.e. engine B tends to identify proteins that other
#' engine don't. But also remember that in your setup each protein needs to
#' be identified by at least two standard engines (only one in open search)


num_peps <- compare_first_sec_W(M$run, "num_peps", "ProteinId", TRUE)
tb <- inner_join(M$data, M$lfq, by = join_by(ProteinId))
num_ids <- tb %>%
  filter(ProteinGroupId != "U") %>%
  select(c(ProteinGroupId, pcoverage_align, ProteinId, num_peps, num_unique_peps)) %>%
  mutate(engine_count = purrr::map_dbl(ProteinGroupId, \(x) {
    x <- split_group_str(x, TRUE, TRUE) %>%
      discard(\(x) x == "U")
    return(length(x))
  }))

# Figure out which engines had the biggest contributions
engine_counts <- num_ids$ProteinGroupId %>%
  lapply(., \(x) split_group_str(x, TRUE)) %>%
  unlist() %>%
  discard(\(x) x == "U") %>%
  table()

#' Assign the values of `vec` into three levels: "low", "medium", "high"
#' based on their locations in the distribution of vec
#'
#' @description
#' low = below 1st quartile
#' medium = betweeen 2nd and 3rd quartiles
#' high = above 3rd
levels_by_quartile <- function(vec) {
  quartiles <- summary(vec)
  case_when(
    vec < quartiles[["1st Qu."]] ~ "low",
    quartiles[["3rd Qu."]] < vec ~ "high",
    .default = "medium"
  )
}

compare_col <- "length_category"
engine_hits <- ta$get_peptide_match_df(
  M$peptide_map_path,
  M$aq_reformat_path
) |> as_tibble()
ENGINES <- names(engine_counts)

engine_hits <- engine_hits |> mutate(length_category = levels_by_quartile(engine_hits$length))
hits <- engine_hits

# Collapse the contingency table
engine_hits
engine_hits <- engine_hits %>%
  select(-PeptideId) %>%
  filter(!is.na(!!compare_col)) %>%
  group_by(!!as.symbol(compare_col)) %>%
  filter(n() > 30) %>%
  ungroup()

# Compute effect size for significant hits, using the odds ratio
# This represents how much more likely the engine is to identify
# peptides from a protein of the given category
# Confidence interval is 95%
length_tests <- chisqNME(engine_hits, ENGINES, compare_col, "engine", "peptide length")
TABLES$engine_length_chi <- length_tests$gt$chi
TABLES$engine_length_contingency <- length_tests$gt$contingency

# An example of what the contingency table for the tests would look like,
# and what the expected values would be under the chi-square test of independence
categories <- unique(engine_hits[[compare_col]])
current <- "comet"
tl <- table(engine_hits[[compare_col]], engine_hits[[current]])
ct <- "low"
show_contigency <- table(engine_hits[[compare_col]] == ct, engine_hits[[current]])
expected <- format_engine_contingency(show_contigency, ct, TRUE)

TABLES$sample_expected <- expected %>%
  gt() %>%
  tab_header(
    title = glue("Engine: {current}, Category: {ct}"),
    subtitle = "Expected values under chi-square in parentheses"
  )


# Correlation between no. identifications by engines and coverage
engine_cor <- cor.test(num_ids$engine_count, num_ids$pcoverage_align)
engine_cor$data.name <- "Correlation between number of identifications by different engines and coverage"
# A weak positive correlation, but statistically significant
# Correlation between peptide number and coverage
n_peps_cor <- cor.test(num_ids$num_unique_peps, num_ids$pcoverage_align)
engine_cor$data.name <- "Correlation between number of identified peptides and coverage"
TABLES$correlation <- gt(bind_rows(htest2tb(engine_cor), htest2tb(n_peps_cor)))
# A weak positive correlation, but statistically significant

leq_intensity1st <- hits$mean_intensity <= quantile(hits$mean_intensity, 0.25)
geq_intensity3nd <- quantile(hits$mean_intensity, 0.75) <= hits$mean_intensity
hits$intensity_class <- map2_chr(leq_intensity1st, geq_intensity3nd, \(x, y) {
  if (x) {
    "low"
  } else if (y) {
    "high"
  } else {
    "medium"
  }
})

intense <- chisqNME(
  tb = hits, var_a_levels = ENGINES,
  var_b_col = "intensity_class", var_a = "engine", var_b = "intensity_class",
  binary = TRUE
)

TABLES$engine_intensity_chi <- intense$gt$chi
TABLES$engine_intensity_contingency <- intense$gt$contingency

save(c(TABLES, GRAPHS), glue("{M$outdir}/engine_category_bias"))
