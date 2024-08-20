options(
  browser = "firefox",
  rlang_backtrace_on_error = "full",
  error = rlang::entrace
)
rlang::global_entrace()

if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

an <- new.env()
reticulate_show_error(
  reticulate::source_python(glue("{args$python_source}/annotate.py"), envir = an)
)
PALETTE <- "ggthemes::excel_Slipstream"
BIG_PALETTE <- "ggthemes::Tableau_20"
TABLES <- list()
GRAPHS <- list()

# Compare with previous results
p_rename <- c(
  ACC = "Accession Number",
  pcoverage_align = "Sequence coverage [%]",
  peptideIds = "Peptide IDs",
  header = "Fasta headers"
)

type_mapping <- list(
  "cytoskeleton component" = "z",
  "pore-forming toxin" = "venom_component",
  "ribosomal protein" = "j",
  "oxido-reductive" = "c",
  "hydrolase" = "g",
  "heat shock proteins" = "o",
  "receptor" = "t",
  "proteases" = "o",
  "transport protein" = "u",
  "non-motor microtubule binding protein" = "z",
  "transferase" = "o",
  "enzymes" = "q",
  "metallopeptidase" = "o",
  "signalling protein" = "t",
  "transcription factors" = "k",
  "transfer/carrier protein" = "u",
  "toxin" = "venom_component",
  "regulatory protein" = "Miscellaneous",
  "metalloendopeptidases" = "o",
  "unknown" = "Function unknown",
  "scaffold/adaptor protein" = "t",
  "photoprotein" = "Miscellaneous",
  "cleavage stimulation factor" = "a",
  "vesicle coat protein" = "u",
  "miscellaneous" = "Miscellaneous",
  "kinase modulator" = "o",
  "small gtpase" = "u",
  "rna splicing factor" = "a",
  "lipoprotein" = "u",
  "g-protein" = "t",
  "translation factors" = "j",
  "cytoplasmic protein" = "Miscellaneous",
  "cytoskeleton component (membrane traffic protein)" = "u",
  "cold shock proteins" = "o",
  "gtpase" = "t",
  "catalyzes" = "Miscellaneous",
  "multiprotein complex" = "Miscellaneous",
  "metalloprotease" = "o",
  "transmembrane protein" = "u",
  "cytoskeleton component (dna metabolism protein)" = "z",
  "transporter" = "u",
  "catalyse" = "Miscellaneous",
  "gtp-binding protein" = "t",
  "motor protein" = "u",
  "chromatin regulator" = "j",
  "calcium-binding protein" = "t",
  "metabolic enzymes" = "Miscellaneous",
  "tight junction" = "w"
)

p_toxins <- read_tsv(glue("{M$wd}/data/reference/previous_toxins.tsv")) %>%
  rename(., any_of(p_rename)) %>%
  select(-contains(" "))

p_ref <- read_tsv(glue("{M$wd}/data/reference/previous_all.tsv")) |>
  rename(any_of(p_rename)) |>
  select(-c(Score, pcoverage_align))

p_all <- read_tsv(glue("{M$wd}/data/reference/previous_all_mq.tsv")) %>%
  rename(., any_of(p_rename)) %>%
  select(-contains(" "), -contains("raw")) |>
  mutate(ACC = map_chr(header, an$id_from_header_str)) |>
  left_join(p_ref, by = join_by(ACC)) |>
  mutate(assigned_COG = map_chr(`Protein type`, \(x) {
    x <- str_to_lower(x)
    find <- lget(type_mapping, x, NULL)
    if (is.null(find)) {
      return("Function unknown")
    }
    lget(M$cog_map, find, find)
  }))


data <- read_tsv(M$data_w_cat_path)

p_grouped <- p_all |>
  group_by(peptideIds) |>
  summarise(
    assigned_COG = modes(assigned_COG, first = TRUE, drop = "Function unknown"),
    pcoverage_align = max(pcoverage_align)
  ) |>
  mutate(source = "previous")

cur_grouped <- data |>
  group_by(GroupUP) |>
  summarise(assigned_COG = modes(assigned_COG, first = TRUE, drop = "Function unknown")) |>
  mutate(assigned_COG = replace(assigned_COG, is.na(assigned_COG), "Function unknown"))

prev_cog_counts <- table(p_grouped$assigned_COG) |>
  table2tb(id_col = "COG") |>
  filter(COG != "Function unknown") |>
  mutate(source = "Previous", prop = n / sum(n))
cog_counts <- table(cur_grouped$assigned_COG) |>
  table2tb(id_col = "COG") |>
  filter(COG != "Function unknown") |>
  mutate(source = "New", prop = n / sum(n))

compare_cogs <- bind_rows(cog_counts, prev_cog_counts) |> simplify_cog(cog_col = "COG")

GRAPHS$compare_cog <- ggplot(compare_cogs, aes(y = prop, x = source, fill = str_wrap(COG, 30))) +
  geom_col() +
  M$default_theme +
  scale_fill_paletteer_d(BIG_PALETTE) +
  ylab("Proportion") +
  xlab("Source") +
  guides(fill = guide_legend("Assigned COG"))

compare_all <- inner_join(p_all, M$data, by = join_by(header)) |>
  mutate(
    pcoverage_align.y = pcoverage_align.y * 100,
    ID = map2_chr(NCBI_ID, UniProtKB_ID, \(x, y) ifelse(is.na(x), y, x))
  )

grouped_cov <- compare_all |>
  group_by(peptideIds.x) |>
  summarize(
    pcoverage_align.x = max(pcoverage_align.x),
    pcoverage_align.y = max(pcoverage_align.y)
  ) |>
  rename(previous = pcoverage_align.x, current = pcoverage_align.y)

cov_longer <- grouped_cov %>%
  mutate(
    diff = current - previous,
  ) %>%
  select(-c(current, previous)) %>%
  pivot_longer(cols = !peptideIds.x) %>%
  mutate(value = round(value, 2))


GRAPHS$shared_cov_bp <- cov_longer %>% ggplot(aes(y = value)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  ylab("Percent change (%)") +
  ggtitle("Increase in shared protein sequence coverage from previous analysis")


compare_long <- compare_all |>
  select(ID, pcoverage_align.x, pcoverage_align.y) |>
  rename(previous = pcoverage_align.x, first = pcoverage_align.y)


GRAPHS$protein_wise_coverage2 <- local({
  compare_long <- grouped_cov |>
    mutate(change = current - previous, Improved = change > 0)
  ggplot(
    compare_long,
    aes(x = previous, y = change, color = Improved)
  ) +
    geom_point() +
    scale_fill_paletteer_d(PALETTE) +
    xlab("Previous coverage (%)") +
    ylab("Change in coverage (%)")
})
GRAPHS$protein_wise_coverage2

p_test <- wilcox.test(grouped_cov$previous, grouped_cov$current, paired = TRUE) |>
  htest2tb(data.name = "Coverage of previous vs current")
if (p_test$p_value < 0.05) {
  p_test <- p_test |> bind_rows(
    wilcox.test(grouped_cov$previous, grouped_cov$current, paired = TRUE, alternative = "less") |>
      htest2tb(data.name = "Previous x current", alternative = "less")
  )
}

dgrouped <- data |>
  filter(!grepl("D", ProteinId)) |>
  group_by(GroupUP) |>
  summarise(pcoverage_align = max(pcoverage_align) * 100) |>
  mutate(source = "current")


p_test <- p_test |> bind_rows(
  wilcox.test(p_grouped$pcoverage_align, dgrouped$pcoverage_align) |>
    htest2tb(data.name = "Previous x current (unpaired)")
)


TABLES$tests <- p_test

save(c(GRAPHS, TABLES), glue("{M$outdir}/comparison_with_previous"))
