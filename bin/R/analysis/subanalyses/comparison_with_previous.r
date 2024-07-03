if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

PALETTE <- "ggthemes::excel_Slipstream"
TABLES <- list()
GRAPHS <- list()

# Compare with previous results
p_rename <- c(
  NCBI_ID = "Accession Number",
  pcoverage_nmatch.prev = "Sequence coverage [%]"
)
p_all <- read_tsv(glue("{M$wd}/data/reference/previous_all.tsv")) %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))
p_toxins <- read_tsv(glue("{M$wd}/data/reference/previous_toxins.tsv")) %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))

has_id <- M$data |>
  filter(!is.na(NCBI_ID)) |>
  mutate(
    NCBI_ID = map_chr(NCBI_ID, \(x) {
      str_replace(x, "\\.1$", "")
    }),
    UniProtKB_ID = map_chr(UniProtKB_ID, \(x) {
      if (is.na(x)) {
        return(x)
      }
      str_replace(x, "\\_.*", "")
    })
  )
uniprot_joined <- inner_join(p_all, has_id, by = join_by(x$NCBI_ID == y$UniProtKB_ID))
compare_all <- inner_join(p_all, has_id, by = join_by(NCBI_ID)) |>
  bind_rows(uniprot_joined) |>
  mutate(
    pcoverage_nmatch = pcoverage_nmatch * 100,
    ID = map2_chr(NCBI_ID, UniProtKB_ID, \(x, y) ifelse(is.na(x), y, x))
  )

cov_longer <- compare_all %>%
  select(contains("coverage"), ID) %>%
  rename(first = pcoverage_nmatch, previous = pcoverage_nmatch.prev) |>
  mutate(
    diff = first - previous,
  ) %>%
  select(-c(first, previous, pcoverage_align)) %>%
  pivot_longer(cols = !ID) %>%
  mutate(value = round(value, 2))

GRAPHS$shared_cov_bp <- cov_longer %>% ggplot(aes(y = value)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  ylab("Percent change (%)") +
  ggtitle("Increase in shared protein sequence coverage from previous analysis")


# For space purposes, drop entries where
# ALL had less than 30% sequence coverage, as well as any entries both second and first passes cannot be compared

compare_long <- compare_all |>
  select(ID, pcoverage_nmatch.prev, pcoverage_nmatch) |>
  rename(previous = pcoverage_nmatch.prev, first = pcoverage_nmatch) |>
  mutate(greater = map2_dbl(first, previous, \(x, y) max(x, y)), previous = -previous) |>
  pivot_longer(cols = -c(ID, greater)) |>
  mutate(is_greater = map2_lgl(value, greater, \(x, y) abs(x) == y))


n_greater <- local({
  previous <- filter(compare_long, name == "previous")$is_greater |> sum()
  first <- filter(compare_long, name == "first")$is_greater |> sum()
  list(previous = previous, first = first)
})


GRAPHS$protein_wise_coverage2 <- local({
  compare_long <- compare_all |>
    select(ID, pcoverage_nmatch.prev, pcoverage_nmatch) |>
    rename(previous = pcoverage_nmatch.prev, first = pcoverage_nmatch) |>
    mutate(change = first - previous, Improved = change < 0)
  ggplot(
    compare_long,
    aes(x = previous, y = change, color = Improved)
  ) +
    geom_point() +
    scale_fill_paletteer_d(PALETTE) +
    xlab("Previous coverage (%)") +
    ylab("Change in coverage (%)")
})



protein_wise_coverage <- local({
  plot <- compare_long |>
    ggplot(aes(x = ID, y = value, fill = name)) +
    geom_bar(stat = "identity", position = "identity") +
    ylab("Coverage (%)") +
    theme(
      axis.text.x = element_blank(),
    ) +
    xlab("Protein") +
    ggtitle("Per-protein sequence coverage") +
    scale_fill_paletteer_d(PALETTE) +
    guides(fill = guide_legend("Run")) +
    scale_y_continuous(
      breaks = seq(-100, 100, by = 20),
      labels = abs(seq(-100, 100, by = 20)),
      limits = c(-100, 100)
    )
  for (i in seq_len(nrow(compare_long))) {
    row <- compare_long[i, ]
    if (row$is_greater) {
      plot <- plot + geom_point(
        x = row$ID, y = row$value, shape = 18,
        show.legend = FALSE, color = "#f14124", alpha = 0.7
      )
    }
  }
  plot
})
GRAPHS$protein_wise_coverage <- protein_wise_coverage

p_test <- wilcox.test(compare_all$pcoverage_nmatch.prev,
  compare_all$pcoverage_nmatch,
  paired = TRUE, alternative = "less"
) |>
  to("data.name", "Coverage of proteins identified in maxquant-only run vs pipeline") |>
  to("alternative", "maxquant-only run was less") |>
  htest2tb()

TABLES$test <- p_test

# Filter out new proteins
compare_toxin <- compare_all %>% filter(NCBI_ID %in% p_toxins$NCBI_ID)
new_proteins <- M$run$first %>% filter(!NCBI_ID %in% compare_all$NCBI_ID)
new_toxins <- new_proteins %>% filter(category == "venom_component")

# TODO: Get the identity of the previous proteins so you
# can
save(c(GRAPHS, TABLES), glue("{M$outdir}/figures/comparison_with_previous"))
