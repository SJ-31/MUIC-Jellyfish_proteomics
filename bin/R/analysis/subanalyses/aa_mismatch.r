if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
library(magrittr)
PALETTE <- "ggthemes::hc_darkunica"
GRAPHS <- list()
TABLES <- list()

# coverage_threshold <- 0.8
data <- M$data |> distinct(ProteinId, .keep_all = TRUE)
nd_data <- get_run("ND_C_indra", glue("{M$wd}/results/ND_C_indra"),
  which = M$chosen_pass
) %>% filter(ProteinId %in% data$ProteinId)
nd_alignments <- get_alignment_data(glue("{M$wd}/results/ND_C_indra"), which = M$chosen_pass)
alignments <- M$alignments

to_keep <- data %>%
  filter(!is.na(MatchedPeptideIds)) %>%
  purrr::pluck("ProteinId")

to_keep_denovo <- data %>%
  filter(!is.na(MatchedPeptideIds)) %>%
  purrr::pluck("MatchedPeptideIds") %>%
  lapply(., str_split_1, pattern = ";") %>%
  unique() %>%
  unlist()

it <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = it)

# Filter out mismatches in the `default` run that arise from engine peptides, in order
# to compare only UPs to the engine peptides of `ND`
alignments$mismatches <- it$only_denovo(alignments$peptides, alignments$mismatches) |> as_tibble()

# Want to find the de novo peptides that match to engine peptides in the same region
alignments$peptides |> group_by()
# Gets 57,544

mismatch_metrics <- map2(
  list(alignments, nd_alignments),
  list(data, nd_data),
  \(x, y) {
    it$classify_mismatches(x$mismatches) %>%
      it$aggregate_mismatches() %>%
      as_tibble() %>%
      filter(ProteinId %in% y$ProteinId)
  }
) %>%
  `names<-`(c("default", "no_denovo"))


merged_default <- inner_join(mismatch_metrics$default, data)
merged <- bind_rows(
  merged_default |> mutate(mode = "default"),
  inner_join(mismatch_metrics$no_denovo, nd_data) %>% mutate(mode = "no_denovo"),
)


default_pep_num <- alignments$peptides |>
  filter(id != ProteinId) |>
  group_by(ProteinId) |>
  summarise(num_unique_peps_d = n())

merged_prot_compare <- inner_join(mismatch_metrics$default,
  mismatch_metrics$no_denovo,
  by = join_by(ProteinId),
  suffix = c("_default", "_no_denovo")
) |>
  filter(ProteinId %in% filter(data, !is.na(MatchedPeptideIds))$ProteinId) |>
  inner_join(select(data, ProteinId, length), by = join_by(ProteinId)) |>
  inner_join(select(nd_data, ProteinId, num_unique_peps), by = join_by(ProteinId)) |>
  rename(num_unique_peps_nd = num_unique_peps) |>
  inner_join(default_pep_num, by = join_by(ProteinId)) |>
  mutate(
    mean_mismatch_d = n_mismatches_default / num_unique_peps_d,
    mean_mismatch_nd = n_mismatches_no_denovo / num_unique_peps_nd
  )

TABLES$n_mismatches_test <- with(
  merged_prot_compare,
  wilcox.test(n_mismatches_default, n_mismatches_no_denovo, "greater", paired = TRUE)
) |>
  htest2tb()

GRAPHS$mismatch_comparison <- compare_vals_x_y(
  "Average mismatches with denovo peptides",
  "Average mismatches without de novo peptides",
  "mean_mismatch_d", "mean_mismatch_nd",
  merged_prot_compare,
  color = "length",
  palette = "viridis::inferno"
)

TABLES$mismatch_metrics <- gt(mismatch_metrics$default)


GRAPHS$conservative_ratio <- ggplot(merged, aes(x = nc_c_ratio, fill = mode)) +
  geom_histogram(position = "identity", alpha = 0.7) +
  xlab("Ratio of non-conservative to conservative mismatches") +
  scale_fill_paletteer_d(PALETTE)


mapped_by <- list()
mapped_by$denovo <- data %>%
  filter(grepl("D", MatchedPeptideIds) | grepl("D", ProteinId))
mapped_by$transcriptome <- data %>%
  filter(grepl("T", MatchedPeptideIds) | grepl("T", ProteinId))

#' How to get metrics for substitutions on de novo peptides directly?
#' Map peptides from the "alignments" file onto the original sequences of the denovo peptides
#' Get the original sequence
denovo_metrics <- reticulate_show_error(it$denovo_mismatch_metrics(
  to_keep,
  to_keep_denovo, M$seq_map_path, M$unmatched_path,
  alignments$peptides, alignments$mismatches
)) %>% lapply(as_tibble)

GRAPHS$replacement_hist <- denovo_metrics$metrics %$%
  gg_numeric_dist(list(
    conservative = n_conservative,
    `non conservative` = n_non_conservative
  ), "hist", position = "identity", alpha = 0.6) +
  xlab("Replacement count") +
  ylab("Frequency") +
  labs(
    title = "Amino acid substitution frequency",
    subtitle = glue("from proteins with > {coverage_threshold} coverage"),
    fill = "Type",
    color = "Type"
  ) + scale_fill_paletteer_d(PALETTE)

TABLES$denovo_metrics <- gt(denovo_metrics$metrics)

test <- wilcox.test(denovo_metrics$metrics$n_mismatches, y = NULL) |>
  to("data.name", "replacement counts") |>
  htest2tb()

TABLES$replacement_htests <- gt(test)


save(c(GRAPHS, TABLES), glue("{M$outdir}/amino_acid_mismatches"))
