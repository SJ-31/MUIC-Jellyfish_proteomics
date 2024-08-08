library("ggVennDiagram")
library("cowplot")
library("gridExtra")
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
GRAPHS <- list()
TABLES <- list()
## -# Coverage metrics
PALETTE <- "ggthemes::Classic_20"
cov_align <- compare_first_sec_L(
  M$run, "pcoverage_align",
  TRUE, "ProteinId"
)
GRAPHS$run_coverage <- pass_density_plot(cov_align, 0.05) + labs(x = "percent coverage")

# Unique proteins to each run
run_uniques <- get_pass_uniques(M$run)
percent_found <- dplyr::bind_cols(
  not_missing(M$run$first),
  not_missing(M$run$sec)
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
tb <- M$data
grouping_metric <- "assigned_COG"
if (grouping_metric %in% colnames(tb)) {
  lfq <- dplyr::select(tb, all_of(grouping_metric), ProteinId) %>% inner_join(., merge_lfq(tb, "mean"))
  apply_over <- tb[[grouping_metric]] %>%
    table() %>%
    discard(., \(x) x < 100) %>%
    names()
  cov_list <- group_list_from_tb(tb,
    v = apply_over, col_from = grouping_metric,
    target_col = "pcoverage_align"
  )
  intensity_list <- group_list_from_tb(lfq, apply_over, grouping_metric, "log_intensity")
  GRAPHS$intensity_categories <- gg_numeric_dist(intensity_list, "boxplot") +
    labs(y = "log intensity", x = grouping_metric) + theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    guides(color = guide_legend(grouping_metric)) + scale_color_paletteer_d(PALETTE)
  GRAPHS$coverage_categories <- gg_numeric_dist(cov_list, "boxplot") +
    labs(y = "coverage (%)", x = grouping_metric) + theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    guides(color = guide_legend(grouping_metric)) + scale_color_paletteer_d(PALETTE)

  with_category <- inner_join(tb, lfq) %>%
    select(ProteinId, log_intensity, !!grouping_metric) %>%
    filter(!is.na(log_intensity)) %>%
    arrange(log_intensity) %>%
    mutate(rank = seq_len(nrow(.)))
  GRAPHS$category_ranks <- with_category %>%
    ggplot(aes(x = rank, y = log_intensity, color = !!as.symbol(grouping_metric))) +
    geom_point() +
    labs(x = "Rank", y = "Log intensity") +
    scale_color_paletteer_d(PALETTE)

  # Top ten most intense proteins
  top_ten <- lfq %>%
    arrange(desc(log_intensity)) %>%
    slice(1:10)
}

# --------------------------------------------------------

## -# Annotation metrics
counts <- list()
counts$first <- get_counts(M$run$first)
counts$sec <- get_counts(M$run$sec)

wanted_cols <- c("num_peps", "num_unique_peps", "pcoverage_align")
# Run paired Wilcox tests between on the metrics defined above, pairing up proteins that
# were identified in both runs
# Do this for all different run parameters to see if results are consistent

venn_helper <- function(lst, palette, caption) {
  ggVennDiagram(lst,
    label_size = 5,
  ) + scale_fill_paletteer_c(palette) +
    labs(caption = caption) +
    theme(
      plot.caption = element_text(face = "bold", size = 20),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 13)
    )
}

cov_compare_helper <- function(tb, palette, caption) {
  ggplot(tb, aes(
    x = pcoverage_align.first,
    y = pcoverage_align.sec, color = length.first
  )) +
    geom_point() +
    paletteer::scale_colour_paletteer_c(palette) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1),
      linetype = 2,
      colour = "black"
    ) +
    ylab("Second pass coverage") +
    xlab("First pass coverage") +
    guides(color = guide_legend("Protein length")) +
    labs(caption = caption) +
    theme(
      plot.caption = element_text(face = "bold", size = 15)
    )
}


peptide_cols <- c("length", "modified")
wilcox <- tibble()
wilcox_peptides <- tibble()
venns <- list()
covs <- list()
pep_lengths <- list()
palettes <- c("ggthemes::Classic Orange", "ggthemes::Classic Blue", "ggthemes::Classic Green", "ggthemes::Classic Red")
dpalettes <- c("PrettyCols::Neon", "PrettyCols::Lucent", "PrettyCols::Light", "PrettyCols::Joyful")
names <- c("default", "msConvert", "Calibrated", "ND")
for (i in seq_along(M$prefixes)) {
  run <- get_run(M$prefixes[[i]], M$all_paths[[i]])
  compare_tb <- inner_join(run$first, run$second, by = join_by(header), suffix = c(".first", ".sec")) |>
    filter(!grepl("D", ProteinId.first))
  w_temp <- pairwise_tests_tb(
    compare_tb, wanted_cols, c("less", "less", "less"),
    \(x, y, ...) wilcox.test(x, y, paired = TRUE, na.rm = TRUE, ...)
  ) |> mutate(param = M$prefixes[[i]])

  first_peps <- get_peptide_data(run$first$peptideIds) |> mutate(pass = "First")
  second_peps <- get_peptide_data(run$second$peptideIds) |> mutate(pass = "Second")
  only_in_sec <- second_peps |> filter(!peptide %in% first_peps$peptide)
  peptides <- bind_rows(first_peps, only_in_sec) |>
    mutate(modified = as.double(modified))
  w_peptides <- lapply(peptide_cols, \(x) {
    name <- glue("first x second {names[i]}, {x}")
    two_sided <- wilcox.test(
      filter(peptides, pass == "First")[[x]],
      filter(peptides, pass == "Second")[[x]]
    ) |> htest2tb(data.name = name)
    one_sided <- wilcox.test(
      filter(peptides, pass == "First")[[x]],
      filter(peptides, pass == "Second")[[x]],
      "less"
    ) |>
      htest2tb(data.name = name, alternative = "first less")
    bind_rows(two_sided, one_sided) |> mutate(param = names[i], var = x)
  })

  pep_lengths[[names[i]]] <- ggplot(peptides, aes(x = length, y = pass, fill = pass)) +
    geom_density_ridges2() +
    scale_fill_paletteer_d(dpalettes[i]) +
    theme(
      axis.title.y = element_blank(), axis.text.y = element_blank(),
      axis.text.x = element_text(size = 13), axis.title.x = element_text(size = 15)
    ) +
    scale_x_continuous(breaks = round(seq(min(peptides$length), max(peptides$length), 5))) +
    labs(caption = names[i]) +
    theme(
      plot.caption = element_text(face = "bold", size = 15)
    )

  wilcox_peptides <- bind_rows(wilcox_peptides, w_peptides)
  wilcox <- bind_rows(wilcox, w_temp)
  header_lst <- list(`1st` = run$first$header, `2nd` = run$second$header)
  venns[[names[i]]] <- venn_helper(header_lst, palettes[i], names[i])
  covs[[names[i]]] <- cov_compare_helper(compare_tb, palettes[i], names[i])
}

GRAPHS$header_overlap <- do.call("plot_grid", c(venns))
if (!interactive()) {
  GRAPHS$pass_cov_comparison <- do.call("grid.arrange", c(covs, ncol = 2))
}

compare_tb <- inner_join(M$run$first, M$run$second, by = join_by(header), suffix = c(".first", ".sec")) |>
  filter(!grepl("D", ProteinId.first))

sub <- substitute_all(
  c("num_peps", "pcoverage_align"),
  c("peptide number", "percent coverage"),
  \(x) gsub("_", " ", x)
)

GRAPHS$peptide_lengths <- do.call("grid.arrange", c(pep_lengths, ncol = 2))

wilcox_peptides <- wilcox_peptides |>
  mutate(
    p_adjust = p.adjust(p_value),
    significant = as.double(p_adjust < 0.05)
  ) |>
  rename(pair = data)

TABLES$wilcox_peptides <- wilcox_peptides |>
  conclude_one_sided() |>
  select(-c(param, var)) |>
  pairwise_conclusion2gt()

TABLES$wilcox_pass_comparison <- wilcox %>%
  rename_with(
    \(vec) {
      rename_helper <- function(x) {
        if (str_detect(x, "sided")) {
          x <- str_replace(x, ".sided", "-sided")
        }
        if (str_detect(x, "p_value")) {
          x <- str_replace(x, "p_value", "p-value")
        }
        str_replace_all(x, "_", " ")
      }
      map_chr(vec, rename_helper)
    },
    everything()
  ) |>
  gt(rowname_col = "metric") %>%
  fmt(., columns = "metric", fns = \(x) map_chr(x, sub)) %>%
  text_case_match(
    "two.sided" ~ "two-sided",
    ".first" ~ "first",
    "_" ~ " pass ",
    .replace = "partial"
  ) %>%
  fmt_number(., columns = contains("p-value"), decimals = 5) %>%
  tab_stubhead(label = "Metric")

# Results per protein
# Evaluate significance of each
per_protein <- tibble(
  metric = rep(wanted_cols, 2),
  type = c(rep("mean", length(wanted_cols)), rep("stdev", length(wanted_cols))),
  first = c(
    get_avg_sd(M$run$first, wanted_cols, \(x) mean(x, na.rm = TRUE)),
    get_avg_sd(M$run$first, wanted_cols, \(x) sd(x, na.rm = TRUE))
  ),
  sec = c(
    get_avg_sd(M$run$sec, wanted_cols, \(x) mean(x, na.rm = TRUE)),
    get_avg_sd(M$run$sec, wanted_cols, \(x) sd(x, na.rm = TRUE))
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
  geom_bar(stat = "identity") +
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())


save(c(GRAPHS, TABLES), glue("{M$outdir}/general_metrics"))
