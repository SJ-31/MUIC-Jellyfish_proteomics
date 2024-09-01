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

record_stats <- function(tb, filter_list, source) {
  stats <- list(source = c(), count = c(), prop = c(), total = c(), criteria = c())
  size <- nrow(tb)
  for (n in names(filter_list)) {
    filter <- filter_list[[n]]
    f <- tb |> filter(filter)
    count <- nrow(f)
    stats$count <- append(stats$count, count)
    stats$source <- append(stats$source, source)
    stats$total <- append(stats$total, size)
    stats$prop <- append(stats$prop, count / size)
    stats$criteria <- append(stats$criteria, n)
  }
  return(as_tibble(stats))
}

# Check if de novo peptides and full length-proteins are ever found together
all_def <- read_tsv(glue("{denovo_dir}/default_prot_all.tsv"))
all_joined <- read_tsv(glue("{denovo_dir}/joined_prot_all.tsv"))
all_joined_pep <- all_joined |> filter(Peptide == Peptide_nd)

def_filters <- with(all_def, list(
  "Matched to both" = n_denovo > 0 & n_full > 0,
  "Matched to de novo only" = n_denovo > 0 & n_full == 0,
  "Matched to DBP only" = n_denovo == 0 & n_full > 0
))

def_stats <- record_stats(all_def, def_filters, "default")

all_joined_filters <- with(all_joined, list(
  "Matched to both" = n_denovo > 0 & n_full > 0,
  "Matched to de novo only" = n_denovo > 0 & n_full == 0,
  "Matched to DBP only" = n_denovo == 0 & n_full > 0
))

clean_peptide_all <- function(x) {
  clean_peptide(x) |>
    str_remove_all("\\.") |>
    str_remove_all("-") |>
    str_remove_all("n")
}

pyt <- reticulate::import("pyteomics.mass")
compare_peps <- all_joined |>
  filter(n_denovo > 0 & n_full == 0) |>
  select(Peptide, Peptide_nd) |>
  mutate(across(everything(), clean_peptide_all),
    n_mismatch = map2_dbl(Peptide, Peptide_nd, \(x, y) {
      stringdist::stringdist(x, y, method = "lv")
    }),
    mass = map_dbl(Peptide, pyt$fast_mass),
    mass_nd = map_dbl(Peptide_nd, pyt$fast_mass),
    mass_diff = abs(mass - mass_nd)
  )

see(compare_peps)

compare_peps$n_mismatch |> hist()

aj_stats <- record_stats(all_joined, all_joined_filters, "ND x default")


same_pep_filters <- with(all_joined_pep, list(
  "Matched to de novo only" = n_denovo > 0 & n_full == 0,
  "Matched to both" = n_denovo > 0 & n_full > 0,
  "Matched to DBP only" = n_denovo == 0 & n_full > 0,
  "Identical set of protein matches" = Proteins == Proteins_nd
))
aj_stats

sp_stats <- record_stats(all_joined_pep, same_pep_filters, "ND x default shared peptides")

to_plot <- bind_rows(def_stats, aj_stats, sp_stats)

TABLES$psm_stats <- to_plot |> gt()

sp_stats

save(c(TABLES, GRAPHS), denovo_dir)
