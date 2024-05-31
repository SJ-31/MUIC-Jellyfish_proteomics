TABLES <- list()
GRAPHS <- list()

# Compare with previous results
p_rename <- c(
  NCBI_ID = "Accession Number",
  pcoverage_nmatch.prev = "Sequence coverage [%]"
)
p_all <- read_tsv(glue("{wd}/data/reference/previous_all.tsv")) %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))
p_toxins <- read_tsv(glue("{wd}/data/reference/previous_toxins.tsv")) %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))

compare_all <- inner_join(p_all, run$first, by = join_by(NCBI_ID)) %>%
  left_join(., run$sec, by = join_by(NCBI_ID), suffix = JOIN_SUFFIX) %>%
  select(c(
    NCBI_ID, header.first, category.first,
    pcoverage_nmatch.prev, pcoverage_nmatch.first, pcoverage_nmatch.sec
  )) %>%
  mutate(
    pcoverage_nmatch.first = pcoverage_nmatch.first * 100,
    pcoverage_nmatch.sec = pcoverage_nmatch.sec * 100
  )

cov_longer <- compare_all %>%
  select(contains("coverage"), NCBI_ID) %>%
  rename_with(., \(x) purrr::map_chr(x, \(y) {
    if (grepl("NCBI_ID", y)) {
      return(y)
    }
    y <- gsub("pcoverage_nmatch", "", y)
    if (y == ".prev") {
      return("previous")
    } else if (y == ".first") {
      return("first")
    }
    return("second")
  })) %>%
  mutate(
    first_diff = first - previous,
    second_diff = second - previous
  ) %>%
  select(-c(first, second, previous)) %>%
  pivot_longer(cols = !NCBI_ID) %>%
  mutate(value = round(value, 2))

GRAPHS$shared_cov_bp <- cov_longer %>% ggplot(aes(y = value, fill = name)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  xlab("Seqence coverage") +
  ylab("Percent change (%)") +
  scale_fill_discrete(name = "Pass", labels = c("first", "second")) +
  ggtitle("Increase in shared protein sequence coverage")


# For space purposes, drop entries where
# ALL had less than 30% sequence coverage, as well as any entries both second and first passes cannot be compared
threshold <- 0
cov_longer <- cov_longer %>%
  group_by(NCBI_ID) %>%
  nest()
# mask <- cov_longer %>% apply(1, \(x) {
#   data <- x$data
#   if (all(data$value < threshold) | any(is.na(data$value))) {
#     return(FALSE)
#   }
#   return(TRUE)
# })
# cov_longer <- cov_longer[mask, ] %>% unnest(cols = c(data))
cov_longer$name <- factor(cov_longer$name, levels = c("first_diff", "second_diff"))

GRAPHS$protein_wise_coverage <- cov_longer %>%
  ggplot(aes(x = NCBI_ID, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Coverage change (%)") +
  xlab(element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  ggtitle("Per-protein sequence coverage change") +
  scale_fill_discrete(name = "Pass", labels = c("first", "second"))

p_test <- wilcox.test(compare_all$pcoverage_nmatch.prev,
  compare_all$pcoverage_nmatch.first,
  paired = TRUE, alternative = "less"
) |>
  to("data.name", "Coverage of proteins identified in maxquant-only run vs pipeline") |>
  to("alternative", "maxquant-only run was less") |>
  htest2Tb()

# Filter out new proteins
compare_toxin <- compare_all %>% filter(NCBI_ID %in% p_toxins$NCBI_ID)
new_proteins <- run$first %>% filter(!NCBI_ID %in% compare_all$NCBI_ID)
new_toxins <- new_proteins %>% filter(category == "venom_component")

# TODO: Get the identity of the previous proteins so you
# can
save(c(GRAPHS, TABLES), glue("{OUTDIR}/figures/comparison_with_previous"))
