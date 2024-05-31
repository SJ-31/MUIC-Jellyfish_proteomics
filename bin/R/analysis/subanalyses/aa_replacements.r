library(magrittr)
GRAPHS <- list()
TABLES <- list()

coverage_threshold <- 0.8
data <- run$first %>% filter(pcoverage_align >= coverage_threshold)
nd_data <- runData("ND_C_indra", glue("{wd}/results/ND_C_indra"),
  which = CHOSEN_PASS
) %>% filter(ProteinId %in% data$ProteinId)
nd_alignments <- alignmentData(glue("{wd}/results/ND_C_indra"), which = CHOSEN_PASS)

unmatched_peptides <- read_tsv(UNMATCHED_PATH) %>%
  mutate(peptideIds = map_chr(peptideIds, \(x) cleanPeptide(x))) %>%
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
reticulate::source_python(glue("{args$python_source}/trace_alignments.py"), envir = it)

replacement_metrics <- map2(
  list(alignments, nd_alignments),
  list(data, nd_data),
  \(x, y) {
    it$classifyReplacements(x$replacements) %>%
      it$aggregateReplacements() %>%
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


GRAPHS$peptides_vs_replacements <- ggplot(merged, aes(
  x = num_unique_peps, y = n_replacements,
  color = mode, alpha = n_conflicts / n_replacements
)) +
  geom_point() +
  xlab("Number of unique peptides") +
  ylab("n replacements") +
  labs(alpha = "n conflicts / n replacements")


test <- cor.test(merged$n_replacements, merged$num_unique_peps) |>
  to("data.name", "number of replacements vs number of unique peptides") |>
  htest2Tb()

TABLES$replacement_metrics <- gt(replacement_metrics$default)

GRAPHS$conservative_ratio <- ggplot(merged, aes(x = nc_c_ratio, fill = mode)) +
  geom_histogram(position = "identity", alpha = 0.7) +
  xlab("Ratio of non-conservative substitutions to conservative")

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
denovo_metrics <- it$denovoReplacementMetrics(
  to_keep,
  to_keep_denovo, SEQ_MAP_PATH,
  alignments$peptides, alignments$replacements
) %>% lapply(as_tibble)


GRAPHS$replacement_hist <- denovo_metrics$metrics %$%
  ggplotNumericDist(list(
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
  )

TABLES$denovo_metrics <- gt(denovo_metrics)

test <- wilcox.test(denovo_metrics$metrics$n_replacements, y = NULL) |>
  to("data.name", "replacement counts")
htest2Tb() |>
  bind_rows(test)

TABLES$replacement_htests <- gt(test)


save(c(GRAPHS, TABLES), glue("{OUTDIR}/figures/amino_acid_replacements"))
