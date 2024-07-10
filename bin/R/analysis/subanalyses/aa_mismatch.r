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

replacement_metrics <- map2(
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

merged_default <- inner_join(replacement_metrics$default, data)
merged <- bind_rows(
  merged_default |> mutate(mode = "default"),
  inner_join(replacement_metrics$no_denovo, nd_data) %>% mutate(mode = "no_denovo"),
)


# TODO: This is kinda weird and misleading, maybe something better...
# GRAPHS$peptides_vs_mismatches <- ggplot(merged, aes(
#   x = num_unique_peps, y = n_mismatches,
#   color = mode, alpha = n_disagreements / n_mismatches
# )) +
#   geom_point() +
#   xlab("Number of unique peptides") +
#   ylab("n mismatches") +
#   labs(alpha = "n disagreements / n mismatches") +
#   scale_color_paletteer_d(PALETTE)

test <- cor.test(merged$n_mismatches, merged$num_unique_peps) |>
  to("data.name", "number of mismatches vs number of unique peptides") |>
  htest2tb()

TABLES$replacement_metrics <- gt(replacement_metrics$default)


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
  htest2tb() |>
  bind_rows(test)

TABLES$replacement_htests <- gt(test)


save(c(GRAPHS, TABLES), glue("{M$outdir}/Figures/amino_acid_mismatches"))
