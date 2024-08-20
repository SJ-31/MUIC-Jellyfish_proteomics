if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
ta <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = ta)

open_search_engines <- c("metamorpheusGPTMD", "msfraggerGPTMD", "msfraggerGlyco")
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

ENGINES <- names(engine_counts)
ENGINES <- replace(ENGINES, ENGINES == "metamorpheusGTPMD", "metamorpheusGPTMD")

all_pepmaps <- lapply(M$all_paths, \(x) {
  c(
    glue("{x}/{M$passes[[1]]}/percolator_peptide_map_all.tsv"),
    glue("{x}/{M$passes[[2]]}/percolator_peptide_map_all.tsv")
  )
}) |> `names<-`(M$params)


all_aq_reformat <- lapply(M$all_paths, \(x) {
  c(
    glue("{x}/{M$passes[[1]]}/Quantify/directlfq.aq_reformat.tsv"),
    glue("{x}/{M$passes[[2]]}/Quantify/directlfq.aq_reformat.tsv")
  )
}) |> `names<-`(M$params)

# Correlation between no. identifications by engines and coverage
engine_cor <- cor.test(num_ids$engine_count, num_ids$pcoverage_align)
engine_cor$data.name <- "Correlation between number of identifications by different engines and coverage"
# A weak positive correlation, but statistically significant
# Correlation between peptide number and coverage
n_peps_cor <- cor.test(num_ids$num_unique_peps, num_ids$pcoverage_align)
engine_cor$data.name <- "Correlation between number of identified peptides and coverage"
TABLES$correlation <- gt(bind_rows(htest2tb(engine_cor), htest2tb(n_peps_cor)))
# A weak positive correlation, but statistically significant

intensity_bias <- list()
matches <- list()
length_bias <- list()
CONTINGENCY <- list()



var_stats <- tibble()
for (i in M$params) {
  pepmap <- all_pepmaps[[i]]
  aq_reformat <- all_aq_reformat[[i]]
  cur_hits <- ta$get_peptide_match_df(
    pepmap[2],
    aq_reformat[2]
  ) |>
    as_tibble()

  cur_hits <- cur_hits |>
    mutate(length_category = levels_by_quartile(cur_hits$length))
  s <- summary_tb(cur_hits$length) |> mutate(param = i, var = "length")
  var_stats <- bind_rows(var_stats, s)


  w_intensity <- cur_hits |> filter(!is.na(mean_intensity))
  leq_intensity1st <- w_intensity$mean_intensity <= quantile(w_intensity$mean_intensity, 0.25)
  geq_intensity3nd <- quantile(w_intensity$mean_intensity, 0.75) <= w_intensity$mean_intensity
  w_intensity$intensity_class <- map2_chr(leq_intensity1st, geq_intensity3nd, \(x, y) {
    if (x) {
      "low"
    } else if (y) {
      "high"
    } else {
      "medium"
    }
  })
  s <- summary_tb(w_intensity$mean_intensity) |> mutate(param = i, var = "PI")
  var_stats <- bind_rows(var_stats, s)
  intense <- chisqNME(
    tb = w_intensity, var_a_levels = ENGINES,
    var_b_col = "intensity_class", var_a = "engine", var_b = "intensity_class",
    binary = TRUE
  )
  # Compute effect size for significant hits, using the odds ratio
  # This represents how much more likely the engine is to identify
  # peptides from a protein of the given category
  # Confidence interval is 95%

  match_type <- chisqNME(
    tb = cur_hits, var_a_levels = ENGINES,
    var_b_col = "match_type", var_a = "engine", var_b = "match_type",
    binary = TRUE
  )

  length_tests <- chisqNME(
    tb = cur_hits, var_a_levels = ENGINES,
    var_b_col = "length_category", var_a = "engine", var_b = "length_category",
    binary = TRUE
  )
  CONTINGENCY[[glue("match_{i}")]] <- match_type$gt$contingency
  CONTINGENCY[[glue("intensity_{i}")]] <- intense$gt$contingency
  CONTINGENCY[[glue("length_{i}")]] <- length_tests$gt$contingency
  matches[[i]] <- match_type$tb$chi |> mutate(param = i)
  intensity_bias[[i]] <- intense$tb$chi |> mutate(param = i)
  length_bias[[i]] <- length_tests$tb$chi |> mutate(param = i)
}

all_intensity_bias <- bind_rows(intensity_bias)
all_length_bias <- bind_rows(length_bias)
all_match_bias <- bind_rows(matches)

or_graph <- function(tb, x, palette) {
  tb |>
    filter(is_significant) |>
    ggplot(aes(x = !!as.symbol(x), y = OR, fill = engine)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~param) +
    M$default_theme +
    geom_errorbar(aes(ymin = OR_lower_ci, ymax = OR_upper_ci),
      position = position_dodge(width = 0.9), width = 0.2
    ) +
    scale_fill_paletteer_d(palette) +
    ylab("OR") +
    annotate("segment",
      x = -Inf, xend = Inf, y = 1, yend = 1, linetype = 2,
      color = "red"
    )
}

GRAPHS$all_intensity_bias <- all_intensity_bias |>
  or_graph("intensity_class", "basetheme::deepblue") + xlab("Intensity class") +
  scale_x_discrete(limits = c("low", "medium", "high"))

GRAPHS$all_length_bias <- all_length_bias |> or_graph(
  "length_category",
  "ggthemr::flat"
) +
  xlab("Length category") + scale_x_discrete(limits = c("low", "medium", "high"))

bias_to_plot <- all_match_bias |>
  filter(param != "ND" & !engine %in% open_search_engines)
GRAPHS$all_match_bias <- or_graph(bias_to_plot, "match_type", "MoMAColors::Warhol") + xlab("Peptide source")

TABLES$match_bias <- all_match_bias |>
  arrange(desc(OR)) |>
  gt()
TABLES$intensity_bias <- all_intensity_bias |>
  arrange(desc(OR)) |>
  gt()
TABLES$length_bias <- all_length_bias |>
  arrange(desc(OR)) |>
  gt()

all_bias <- bind_rows(
  mutate(all_match_bias, var = "source"),
  mutate(all_intensity_bias, var = "intensity"),
  mutate(all_length_bias, var = "length")
) |> mutate(category = dplyr::coalesce(length_category, match_type, intensity_class))

grouped_bias <- all_bias |>
  group_by(engine, var, category) |>
  nest() |>
  mutate(
    mean_OR = map_dbl(data, \(x) mean(x$OR, na.rm = TRUE)),
    sd_OR = map_dbl(data, \(x) sd(x$OR, na.rm = TRUE))
  )
summary(grouped_bias$sd_OR)

# TODO: Can you do this with logistic regression


attr(GRAPHS$all_intensity_bias, "width") <- 18
attr(GRAPHS$all_length_bias, "width") <- 18
attr(GRAPHS$all_match_bias, "width") <- 18


TABLES$bias_stats <- var_stats |>
  group_by(var) |>
  summarize(across(is.double, mean)) |>
  gt() |>
  cols_label(sd = "Standard deviation", var = "Variable") |>
  cols_label_with(fn = str_to_title)

save(c(TABLES, GRAPHS), glue("{M$outdir}/engine_category_bias"))
save(CONTINGENCY, glue("{M$outdir}/engine_category_bias/contingency"))
