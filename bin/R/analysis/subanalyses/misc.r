library("paletteer")
library("ggVennDiagram")
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
TABLES <- list()
GRAPHS <- list()
# --------------------------------------------------------
# Investigating trends in missing quantification
current <- M$data |> inner_join(M$lfq, by = join_by(ProteinId))
cq <- c("directlfq", "maxlfq")
missing_quant_tests <- list()
for (q in cq) {
  noq <- current %>% filter(is.na(!!as.symbol(glue("{cq}_mean"))))
  hasq <- current %>% filter(!ProteinId %in% noq$ProteinId)
  missing_quant_tests[[q]] <- wilcox.test(noq$num_peps, hasq$num_peps,
    alternative = "l"
  )
}
capture.output(missing_quant_tests, file = glue("{M$outdir}/missing_quantification_tests.txt"))
rm(noq)
rm(hasq)
# --------------------------------------------------------

fragment_names <- c("fragment", "partial")
fragment_regex <- paste0(fragment_names, collapse = "|")

unique_entries <- flatten_by(M$data$entry_name, ";") |>
  unique() |>
  discard(\(x) str_detect(x, "-DENOVO|-TRANSCRIPTOME"))

fragments <- M$data |>
  mutate(entry_name = str_to_lower(entry_name)) |>
  filter(grepl(fragment_regex, entry_name))

grouped <- read_tsv(M$data_w_cat_path) |>
  group_by(GroupUP) |>
  summarize(
    assigned_COG = paste0(unique(assigned_COG), collapse = ";"),
    entry_name = paste0(entry_name, collapse = ";"),
    size = n()
  ) |>
  mutate(entry_name = map_chr(entry_name, split_unique_join)) |>
  arrange(desc(size))
TABLES$grouped_cog_sizes <- grouped



# Get Count of fragments


# Just to check if group assignments are correct
cog_evidence <- read_tsv(glue("{M$chosen_path}/Analysis/cog_assignment_evidence.tsv"))
venom_evidence <- cog_evidence |> filter(Group == "venom_component")
venom_gos <- inner_join(venom_evidence, read_tsv(M$go_reference), by = join_by(x$Evidence == y$GO_IDs))

# ----------------------------------------
# Merging results
nd_merged_path <- glue("{M$wd}/results/ND_MERGED")
nd_run <- get_run(M$prefixes[[4]], M$ndpath)
ndm_run <- get_run(M$prefixes[[4]], nd_merged_path)


all_tests <- tibble()
stats <- tibble()
joined <- lapply(c("first", "second"), \(x) {
  tmp_joined <- inner_join(nd_run[[x]], ndm_run[[x]], by = join_by(header), suffix = c(".nd", ".ndm")) |>
    mutate(pass = x)
  test <- with(tmp_joined, wilcox.test(pcoverage_align.nd,
    pcoverage_align.ndm,
    paired = TRUE, alternative = "less"
  )) |>
    htest2tb(data.name = "ND x ND merged", alternative = "ND less") |>
    mutate(pass = x)
  all_tests <<- bind_rows(all_tests, test)
  tmp_stats <- with(tmp_joined, tibble(
    equals = sum(pcoverage_align.nd == pcoverage_align.ndm),
    nd_greater = sum(pcoverage_align.nd > pcoverage_align.ndm),
    ndm_greater = sum(pcoverage_align.nd < pcoverage_align.ndm),
    pass = x
  ))
  stats <<- bind_rows(stats, tmp_stats)

  tmp_joined
}) |>
  bind_rows()

stats <- mutate(stats, prop_ndm_greater = ndm_greater / (ndm_greater + equals + nd_greater)) |> gt()

TABLES$nd_ndm_comparison <- stats

GRAPHS$nd_ndm_comparison <- compare_vals_x_y(
  "ND", "ND merged", "pcoverage_align.nd",
  "pcoverage_align.ndm", joined,
  color = "pass", continuous = FALSE, segment_color = "black"
) +
  M$default_theme
attr(GRAPHS$nd_ndm_comparison, "width") <- 15

# ----------------------------------------

save(c(TABLES, GRAPHS), glue("{M$outdir}/misc"))
# ----------------------------------------
# De novo
denovo_dir <- glue("{M$outdir}/denovo_matches")
all_d <- read_tsv(glue("{denovo_dir}/COMPLETE_final.tsv"))
all_t <- read_tsv(glue("{denovo_dir}/denovo_ND_hits-TRANSCRIPTOME.tsv"))


choose_best <- function(tb) {
  group_by(tb, query) |>
    nest() |>
    mutate(data = lapply(data, \(x) {
      arrange(x, desc(similarity)) |> slice(1)
    })) |>
    unnest(cols = data) |>
    ungroup()
}

from_engines <- all_d |>
  filter(from_engine_ids) |>
  choose_best()

from_ds <- all_d |>
  filter(!from_engine_ids) |>
  choose_best()

matched_denovo_ids <- {
  run <- get_run("C_indra", M$path)
  unique(
    flatten_by(run$second$MatchedPeptideIds, ";")
  ) |> discard(\(x) str_detect(x, "T"))
}

combined <- bind_rows(from_engines, filter(from_ds, !query %in% from_engines$query)) |>
  mutate(is_in_matched = ProteinId %in% matched_denovo_ids)
# Those that aren't matched

D_GRAPHS <- list()

D_GRAPHS$denovo_match_hist <- combined |> ggplot(aes(x = similarity, fill = is_in_matched)) +
  geom_histogram(binwidth = 0.02) +
  M$default_theme +
  ylab("Count") +
  xlab("Levenshtein similarity") +
  guides(fill = guide_legend("Was matched to DBP")) +
  scale_fill_paletteer_d("MoMAColors::Koons")

D_GRAPHS$denovo_match_hist

high_score <- (combined$similarity > 0.7) |> sum()
ratio <- high_score / nrow(combined)

D_GRAPHS$denovo_match_hist_t <- all_t |> ggplot(aes(x = similarity)) +
  geom_histogram(binwidth = 0.02) +
  M$default_theme +
  ylab("Count") +
  xlab("Levenshtein similarity") +
  guides(fill = guide_legend("From engine peptides")) +
  scale_fill_paletteer_d("MoMAColors::Koons")



save(D_GRAPHS, denovo_dir)
