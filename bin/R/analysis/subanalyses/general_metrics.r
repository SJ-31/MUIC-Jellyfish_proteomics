library("ggVennDiagram")
library("cowplot")
library("ggpattern")
library("ggbeeswarm")
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
GRAPHS$run_coverage <- pass_density_plot(cov_align, 0.05) + labs(x = "percent coverage") + M$default_theme

TABLES$run_stats <- get_run_stats(read_tsv(M$data_path)) |>
  gt() |>
  fmt_number(decimals = 0) |>
  cols_label(value = "Value")


stat_cols <- c("ID_method", "inferred_by", "source")
id_stats <- read_tsv(M$data_w_cat_path) |>
  mutate(source = case_when(
    str_detect(ProteinId, "P") ~ "DBP",
    str_detect(ProteinId, "T") ~ "transcriptome",
    .default = "denovo"
  )) |>
  group_by(GroupUP) |>
  summarise(across(all_of(stat_cols), \(x) modes(x, first = TRUE)))

GRAPHS$id_method_source <- id_stats |> ggplot(aes(x = ID_method, fill = source)) +
  geom_bar(position = "dodge") +
  geom_text(aes(label = after_stat(count)),
    stat = "count",
    position = position_dodge(width = 0.9), vjust = 0.1,
    fontface = "bold", size = 5
  ) +
  M$default_theme +
  guides(fill = guide_legend("Protein group type")) +
  ylab("Count") +
  xlab("Identification method") +
  scale_fill_paletteer_d("ltc::fernande")


percent_found <- dplyr::bind_cols(
  not_missing(M$run$first),
  not_missing(M$run$sec)
) |>
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
  scale_fill_discrete("Pass") +
  M$default_theme

# Check if coverage and intensity differs significantly between protein groups
# for confirmation only (we expect them to differ)
grouping_metric <- "assigned_COG"
if (file.exists(M$data_w_cat_path)) {
  w_cat <- read_tsv(M$data_w_cat_path)
  lfq <- merge_lfq(inner_join(w_cat, M$lfq, by = join_by(ProteinId)), "mean")
  tb <- w_cat %>%
    inner_join(lfq, by = join_by(ProteinId)) |>
    simplify_cog() |>
    mutate(assigned_COG = replace(assigned_COG, is.na(assigned_COG), "Function unknown")) |>
    group_by(GroupUP) |>
    summarise(
      assigned_COG = modes(assigned_COG),
      size = n(),
      across(where(is.double), \(x) mean(x, na.rm = TRUE))
    )
  apply_over <- tb[[grouping_metric]] %>%
    table() %>%
    discard(., \(x) x < 100) %>%
    names()
  cov_list <- group_list_from_tb(tb,
    v = apply_over, col_from = grouping_metric,
    target_col = "pcoverage_align"
  )
  intensity_list <- group_list_from_tb(tb, apply_over, grouping_metric, "log_intensity")
  GRAPHS$intensity_categories <- gg_numeric_dist(intensity_list, "boxplot") +
    labs(y = "log intensity", x = grouping_metric) + theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    guides(color = guide_legend(grouping_metric)) + scale_color_paletteer_d(PALETTE) +
    M$default_theme
  GRAPHS$coverage_categories <- gg_numeric_dist(cov_list, "boxplot") +
    labs(y = "Coverage (%)", x = grouping_metric) +
    M$default_theme +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_text(face = "bold")
    ) +
    guides(fill = guide_legend("Assigned COG")) +
    scale_fill_paletteer_d(PALETTE)
  attr(GRAPHS$coverage_categories, "width") <- 18

  with_category <- tb |>
    select(GroupUP, log_intensity, !!grouping_metric, size) %>%
    filter(!is.na(log_intensity) & !is.na(!!grouping_metric)) %>%
    arrange(log_intensity) %>%
    mutate(rank = seq_len(nrow(.)))

  with_category$assigned_COG <- factor(with_category$assigned_COG)
  category_ranks <- with_category %>%
    ggplot(aes(x = rank, y = log_intensity, color = !!as.symbol(grouping_metric))) +
    geom_point() +
    labs(x = "Rank", y = "Log intensity") +
    scale_color_paletteer_d(PALETTE) +
    M$default_theme +
    guides(color = guide_legend("Assigned COG"))

  category_box <- with_category |>
    filter(assigned_COG != "Function unknown") |>
    ggplot(aes(y = log_intensity, fill = str_wrap(assigned_COG, 30))) +
    geom_boxplot() +
    M$default_theme +
    scale_fill_paletteer_d(PALETTE) +
    ylab("Log intensity") +
    theme(axis.text.x = element_blank()) +
    guides(fill = guide_legend("Assigned COG"))

  cog_hist <- tb |>
    filter(assigned_COG != "Function unknown") |>
    ggplot(aes(x = size, fill = assigned_COG)) +
    geom_histogram() +
    scale_y_log10() +
    M$default_theme +
    scale_fill_paletteer_d(PALETTE) +
    guides(fill = "none") +
    ylab("Log 10 count") +
    xlab("Protein group size")

  GRAPHS$cog_hist_intensity <- cowplot::plot_grid(cog_hist, category_box, rel_widths = c(1, 1), labels = "AUTO")
  # Top ten most intense proteins
  attr(GRAPHS$cog_hist_intensity, "width") <- 20
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
  compare_vals_x_y("First pass coverage",
    "Second pass coverage",
    "pcoverage_align.first",
    "pcoverage_align.sec",
    tb,
    color = "length.first",
    palette = palette,
  ) + labs(caption = caption) +
    theme(
      plot.caption = element_text(face = "bold", size = 15)
    ) +
    guides(color = guide_legend("Protein length"))
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

covs$msConvert <- covs$msConvert + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(), axis.title.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), axis.title.y = element_blank()
)
covs$default <- covs$default + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(), axis.title.x = element_blank()
)
covs$ND <- covs$ND + theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(), axis.title.y = element_blank()
)
GRAPHS$header_overlap <- do.call("plot_grid", c(venns)) + M$default_theme
if (!interactive()) {
  GRAPHS$pass_cov_comparison <- do.call("grid.arrange", c(covs, ncol = 2))
  attr(GRAPHS$pass_cov_comparison, "width") <- 16
}

compare_tb <- inner_join(M$run$first, M$run$second, by = join_by(header), suffix = c(".first", ".sec")) |>
  filter(!grepl("D", ProteinId.first))

sub <- substitute_all(
  c("num_peps", "pcoverage_align"),
  c("peptide number", "percent coverage"),
  \(x) gsub("_", " ", x)
)

GRAPHS$peptide_lengths <- do.call("grid.arrange", c(pep_lengths, ncol = 2)) + M$default_theme

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
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank()) +
  M$default_theme

# Format files for easy reading
wanted_read <- c(
  "header", "organism", "lineage", "ID_method", "inferred_by",
  "PANTHER", "interpro_accession", "interpro_description", "interpro_pathways",
  "interpro_db", "PFAMs", "eggNOG_description",
  "entry_name",
  "length", "mass", "num_peps", "num_unique_peps", "GO_counts", "pcoverage_align", "assigned_COG", "GO_IDs"
)

an_dir <- glue("{M$chosen_path}/Analysis")

wcat <- read_tsv(M$data_w_cat_path) |>
  mutate(
    type = case_when(
      str_detect(ProteinId, "D") | str_detect(MatchedPeptideIds, "D") ~ "De novo peptide",
      str_detect(ProteinId, "T") | str_detect(MatchedPeptideIds, "T") ~ "Transcriptome",
      .default = "DBP, standard search"
    ),
    folded_match = case_when(str_detect(ProteinId, "P") ~ "Matched with DBP", .default = "Standalone")
    # Folded match is true if the
    # actual protein is a DBP, and the UP was matched to them
  )

cols <- c("type", "folded_match", "assigned_COG", "ID_method")
wcat_g <- wcat |>
  group_by(GroupUP) |>
  summarize(across(all_of(cols), \(x) modes(x, first = TRUE)))

special <- local({
  m <- filter(wcat_g, ID_method != "standard") |>
    mutate(
      type = "Open search", folded_match =
        case_when(ID_method == "both" ~ "Matched with DBP", .default = "Standalone")
    )
  t <- wcat_g |> filter(type %in% c("De novo peptide", "Transcriptome"))
  bind_rows(m, t)
}) |>
  simplify_cog()

to_plot_dbp <- local({
  tb <- wcat_g |>
    filter(ID_method == "standard" & type == "DBP, standard search") |>
    simplify_cog()
  n <- nrow(tb)
  tb |>
    group_by(assigned_COG) |>
    summarise(prop = n() / n) |>
    mutate(type = "DBP, standard search", folded_match = "DBP")
})

to_plot <- special |>
  group_by(type, folded_match) |>
  select(type, folded_match, assigned_COG) |>
  mutate(size = n()) |>
  group_by(type, folded_match, assigned_COG) |>
  summarise(size_cog = n(), size = unique(size)) |>
  ungroup() |>
  mutate(prop = size_cog / size) |>
  ungroup() |>
  bind_rows(to_plot_dbp)

cog_helper <- function(tb) {
  ggplot(tb, aes(x = type, fill = str_wrap(assigned_COG, 30), y = prop)) +
    geom_col() +
    M$default_theme +
    scale_fill_paletteer_d(PALETTE) +
    guides(fill = guide_legend("Assigned COG"))
}


GRAPHS$type_plot <- to_plot |>
  cog_helper() +
  guides(fill = guide_legend("Assigned COG")) +
  facet_wrap(~folded_match, scales = "free") +
  xlab("Protein group type") +
  ylab("Proportion")
# scale_x_discrete(limits = c("Open search", "Transcriptome", "De novo peptide"))

attr(GRAPHS$type_plot, "width") <- 18

grouped <- wcat |>
  group_by(GroupUP) |>
  summarize(across(all_of(wanted_read), \(x) paste0(x, collapse = ";"))) |>
  mutate(across(all_of(wanted_read), \(x) map_chr(x, split_unique_join))) |>
  left_join(select(tb, GroupUP, log_intensity), by = join_by(GroupUP)) |>
  relocate(c(entry_name, pcoverage_align, log_intensity), .before = everything())


grouped |>
  write_tsv(glue("{an_dir}/grouped_annotations_only.tsv"))

wcat |>
  select(all_of(wanted_read)) |>
  write_tsv(glue("{an_dir}/annotations_only.tsv"))




save(c(GRAPHS, TABLES), glue("{M$outdir}/general_metrics"))
