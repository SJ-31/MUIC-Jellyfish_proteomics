if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
library(magrittr)
PALETTE <- "ggthemes::hc_darkunica"
GRAPHS <- list()
TABLES <- list()

coverage_threshold <- 0.8
data <- M$data %>% filter(pcoverage_align >= coverage_threshold)
nd_data <- get_run("ND_C_indra", glue("{M$wd}/results/ND_C_indra"),
  which = M$chosen_pass
) %>% filter(ProteinId %in% data$ProteinId)
nd_alignments <- get_alignment_data(glue("{M$wd}/results/ND_C_indra"), which = M$chosen_pass)
alignments <- M$alignments

unmatched_peptides <- read_tsv(M$unmatched_path) %>%
  mutate(peptideIds = map_chr(peptideIds, \(x) clean_peptide(x))) %>%
  select(ProteinId, peptideIds)

to_keep <- data %>%
  filter(!is.na(MatchedPeptideIds)) %>%
  purrr::pluck("ProteinId")

to_keep_denovo <- data %>%
  filter(!is.na(MatchedPeptideIds)) %>%
  purrr::pluck("MatchedPeptideIds") %>%
  lapply(., str_split_1, pattern = ";") %>%
  unique() %>%
  unlist() %>%
  c(., unmatched_peptides$ProteinId)

it <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = it)

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

merged_prot_compare <- inner_join(mismatch_metrics$default,
  mismatch_metrics$no_denovo,
  by = join_by(ProteinId),
  suffix = c("_default", "_no_denovo")
) |>
  filter(ProteinId %in% filter(data, !is.na(MatchedPeptideIds))$ProteinId) |>
  inner_join(select(data, ProteinId, length), by = join_by(ProteinId)) |>
  inner_join(select(nd_data, ProteinId, num_unique_peps), by = join_by(ProteinId)) |>
  rename(num_unique_peps_nd = num_unique_peps) |>
  inner_join(select(data, ProteinId, num_unique_peps), by = join_by(ProteinId)) |>
  rename(num_unique_peps_d = num_unique_peps) |>
  mutate(
    mismatch_pr_default = n_mismatches_default / num_unique_peps_d,
    mismatch_pr_no_denovo = n_mismatches_no_denovo / num_unique_peps_nd
  )

# TODO: Because there are so many denovo peptides matched, this doesn't show what you want
GRAPHS$mismatch_comparison <- ggplot(
  merged_prot_compare,
  aes(x = log(mismatch_pr_default), y = log(mismatch_pr_no_denovo), color = length)
) +
  geom_point() +
  paletteer::scale_colour_paletteer_c("viridis::inferno") +
  ylab("Mismatch frequency without de novo peptides") +
  xlab("Mismatch frequency with de novo peptides") +
  geom_segment(aes(x = 0, y = 0, xend = 4, yend = 4))




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
mapped_by$unknown <- data %>%
  filter(grepl("U", MatchedPeptideIds) | grepl("T", ProteinId) | Group == "U")

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
