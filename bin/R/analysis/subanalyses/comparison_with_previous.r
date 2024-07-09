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
  pcoverage_align.prev = "Sequence coverage [%]"
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
    pcoverage_align = pcoverage_align * 100,
    ID = map2_chr(NCBI_ID, UniProtKB_ID, \(x, y) ifelse(is.na(x), y, x))
  )

cov_longer <- compare_all %>%
  select(contains("coverage"), ID) %>%
  rename(first = pcoverage_align, previous = pcoverage_align.prev) |>
  mutate(
    diff = first - previous,
  ) %>%
  select(-c(first, previous, pcoverage_nmatch)) %>%
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
  select(ID, pcoverage_align.prev, pcoverage_align) |>
  rename(previous = pcoverage_align.prev, first = pcoverage_align)

GRAPHS$protein_wise_coverage2 <- local({
  compare_long <- compare_all |>
    select(ID, pcoverage_align.prev, pcoverage_align) |>
    rename(previous = pcoverage_align.prev, first = pcoverage_align) |>
    mutate(change = first - previous, Improved = change > 0)
  ggplot(
    compare_long,
    aes(x = previous, y = change, color = Improved)
  ) +
    geom_point() +
    scale_fill_paletteer_d(PALETTE) +
    xlab("Previous coverage (%)") +
    ylab("Change in coverage (%)")
})

p_test <- wilcox.test(compare_all$pcoverage_align.prev,
  compare_all$pcoverage_align,
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
