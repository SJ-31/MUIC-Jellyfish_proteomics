if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
TABLES <- list()
GRAPHS <- list()

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


GRAPHS$denovo_match_hist <- combined |> ggplot(aes(x = similarity, fill = is_in_matched)) +
  geom_histogram(binwidth = 0.02) +
  M$default_theme +
  ylab("Count") +
  xlab("Levenshtein similarity") +
  guides(fill = guide_legend("Was matched to DBP")) +
  scale_fill_paletteer_d("MoMAColors::Koons")

GRAPHS$denovo_match_hist

high_score <- (combined$similarity > 0.7) |> sum()
ratio <- high_score / nrow(combined)

GRAPHS$denovo_match_hist_t <- all_t |> ggplot(aes(x = similarity)) +
  geom_histogram(binwidth = 0.02) +
  M$default_theme +
  ylab("Count") +
  xlab("Levenshtein similarity") +
  guides(fill = guide_legend("From engine peptides")) +
  scale_fill_paletteer_d("MoMAColors::Koons")

engine_files <- list(
  identipy = "Identipy/identipy_all_pins.temp",
  comet = "Comet/comet_all_pins.temp",
  msfragger = "MsFragger/fragger_all_pins.temp",
  msgf = "MSGF",
  metamorpheus = "Metamorpheus/metamorpheus_AllPSMs_FormattedForPercolator.tab",
  tide = "Tide/tide_search.target.txt"
)

all_def <- read_tsv(glue("{denovo_dir}/default_prot_all.tsv")) |>
  mutate(
    type = case_when(
      n_denovo > 0 & n_full > 0 ~ "Matched to both",
      n_denovo > 0 & n_full == 0 ~ "Matched to de novo only",
      n_denovo == 0 & n_full > 0 ~ "Matched to DBP only"
    ),
    data = "default"
  )
all_joined <- read_tsv(glue("{denovo_dir}/joined_prot_all.tsv"))


# Check if de novo peptides and full length-proteins are ever found together
together <- all_def |> filter(n_full > 0 & n_denovo > 0)
(together$engine |> table()) / nrow(together)
together$n_denovo |> mean()


# Proportion of unique shared spectra where de novo peptides were matched in default, but no DBP was matched in default
matches <- all_joined |> filter(n_denovo > 0)
no_dbp_matches <- all_joined |> filter(n_denovo > 0 & n_full == 0)
no_dbp_matches$n_denovo |> mean()

same_peptides <- all_joined |>
  filter(Peptide == Peptide_nd) |>
  mutate(
    type = case_when(
      n_denovo > 0 ~ "De novo matched in default",
      Proteins == Proteins ~ "Same DBP match",
      .default = "Different DBP match"
    ),
    data = "ND & default joined (same peptides)"
  )


all_joined <- mutate(all_joined,
  type =
    case_when(
      n_denovo > 0 & n_full == 0 ~ "Only de novo matched in default",
      n_denovo > 0 & n_full > 0 ~ "De novo and DBP matched in default",
      n_denovo == 0 ~ "DBP matched only"
    ),
  data = "ND & default joined"
)

to_plot <- bind_rows(all_joined, all_def, same_peptides)


ggplot(to_plot, aes(x = data, fill = type)) +
  geom_bar() +
  scale_y_log10()

to_plot$type |> table()


TABLES$psm_stats <- table(to_plot$type, to_plot$data) |>
  as.data.frame() |>
  as_tibble() |>
  relocate(Var2, .before = everything()) |>
  group_by(Var2) |>
  mutate(proportion = round(Freq / sum(Freq), 2), ) |>
  filter(Freq > 0) |>
  ungroup() |>
  mutate(Index = seq_len(n())) |>
  relocate(Index, .before = everything()) |>
  gt() |>
  cols_label(
    Freq = "Count", proportion = "Proportion",
    Var2 = "PSM source", Var1 = "Type"
  )

same_pep_denovo <- all_joined |> filter(Peptide == Peptide_nd & n_denovo > 0)
same_pep <- all_joined |> filter(Peptide == Peptide_nd & Proteins == Proteins_nd)

TABLES$denovo_and_full <- glue("
Number of PSMs where de novo peptides and full DBPs were matched together: {nrow(together)}
Total: {nrow(all_def)}
Proportion of above out of all PSMs: {nrow(together)/nrow(all_def)}

With respect to ND and default
total joined: {nrow(all_joined)}
1) n unique shared spectra where de novo peptides were matched: {nrow(matches)}, prop: {nrow(matches)/nrow(all_joined)}
2) n unique shared spectra where de novo peptides were matched, but no DBPs matched
in default: {nrow(no_dbp_matches)}, prop: {nrow(no_dbp_matches)/nrow(all_joined)}
Prop of 2) in 1): {nrow(no_dbp_matches)/nrow(matches)}

")


save(c(TABLES, GRAPHS), denovo_dir)
