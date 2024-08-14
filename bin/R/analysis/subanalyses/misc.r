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

# Show that property of every protein in a Percolator group being matched
# by the exact same peptides gets lost when creating new groups via union-find
# data <- M$data %>% inner_join(M$taxa_tb, by = join_by(ProteinId))
# nested <- data |>
#   group_by(Group) |>
#   nest() |>
#   mutate(
#     unique_peptides = map_dbl(data, \(x) {
#       x$peptideIds |>
#         unique() |>
#         length()
#     }),
#     size = map_dbl(data, \(x) nrow(x))
#   ) |>
#   arrange(desc(size))

# ----------------------------------------
# Identify enriched terms by intensity
filter_intensity <- function(tb, quantile = ? Character()) {
  if (quantile == "first") {
    filterFun <- \(x) filter(x, x$log_intensity <= quantile(x$log_intensity, 0.25))
  } else if (quantile == "second") {
    filterFun <- \(x) {
      filter(x, x$log_intensity > quantile(x$log_intensity, 0.25) &
        x$log_intensity < quantile(x$log_intensity, 0.75))
    }
  } else if (quantile == "third") {
    filterFun <- \(x) filter(x, x$log_intensity >= quantile(x$log_intensity, 0.75))
  }
  tb |>
    filterFun() |>
    purrr::pluck("ProteinId")
}

data <- M$data |>
  inner_join(M$lfq, by = join_by(ProteinId)) |>
  distinct(ProteinId, .keep_all = TRUE)
if (!file.exists(glue("{M$ontologizer_path}/high_intensity.tsv"))) {
  ont <- new.env()
  reticulate::source_python(glue("{M$python_source}/ontologizer_wrapper.py"), envir = ont)
  by_intensity <- merge_lfq(data, "mean") %>%
    inner_join(., dplyr::select(data, ProteinId, GO_IDs), by = join_by(ProteinId)) |>
    filter(!is.na(log_intensity))
  O <- ont$Ontologizer(by_intensity, M$ontologizer_exec, M$go_path)
  groups <- list(
    low_intensity = filter_intensity(by_intensity, "first"),
    medium_intensity = filter_intensity(by_intensity, "second"),
    high_intensity = filter_intensity(by_intensity, "third")
  )
  params <- list(`-m` = "Bonferroni-Holm")
  enriched_intensity <- O$runAll(groups, params)

  lmap(enriched_intensity, \(x) {
    write_tsv(x[[1]], glue("{M$ontologizer_path}/{names(x)}.tsv"))
  })
  plot <- ggplot(by_intensity, aes(x = log_intensity)) +
    geom_histogram(fill = "#69d2e7") +
    xlab("log 10 intensity")
  ggsave(glue("{M$outdir}/intensity_histogram.svg"), plot)
} else {
  intensity <- lapply(c("low", "medium", "high"), \(x) read_tsv(glue("{M$ontologizer_path}/{x}_intensity.tsv")))
}

# ----------------------------------------
# Enrich terms based on modifications
mod_names <- c("Met_ox", "Nterm_acetyl", "Lys_acetyl")
if (!file.exists(glue("{M$ontologizer_path}/met_ox.tsv"))) {
  has_mods <- read_tsv(M$percolator_all_path) |>
    filter(ProteinId %in% data$ProteinId) |>
    filter(!is.na(mods)) |>
    fix_all_mods()

  all_mods <- flatten_by(has_mods$mods, ";")
  mod_table <- all_mods |>
    map_chr(\(x) str_extract(x, "(.*\\|.*)\\|[1-9]+", group = 1)) |>
    discard(\(x) str_detect(x, "229.16")) |> # This shouldn't be here...
    table()

  has_mods_ids <- lapply(c("Met\\|15.99", "Nterm\\|42.0", "Lys\\|42.01"), \(x) {
    filter(has_mods, grepl(x, mods)) |> pluck("ProteinId")
  }) |>
    `names<-`(mod_names)
  ont <- new.env()

  reticulate::source_python(glue("{M$python_source}/ontologizer_wrapper.py"), envir = ont)
  O <- ont$Ontologizer(data, M$ontologizer_exec, M$go_path)
  params <- list(`-m` = "Bonferroni-Holm")

  enriched_mods <- O$runAll(has_mods_ids, params)
  lmap(enriched_mods, \(x) {
    write_tsv(x[[1]], glue("{M$ontologizer_path}/{names(x)}.tsv"))
  })
  enriched_mods <- lapply(enriched_mods, as_tibble)
} else {
  enriched_mods <- lapply(mod_names, \(x) {
    read_tsv(glue("{M$ontologizer_path}/{x}.tsv")) |>
      mutate(subset = str_replace(x, "_", " "))
  })
}

# ----------------------------------------
# Analyses for intensity-enriched terms
intensities <- c("high", "medium", "low")
intensity_tbs <- lapply(
  intensities,
  \(x) {
    read_tsv(glue("{M$ontologizer_path}/{x}_intensity.tsv")) |>
      mutate(subset = glue("{x} intensity"))
  }
) |>
  `names<-`(intensities)

intensity_vecs <- intensity_tbs |>
  lapply(\(x) {
    x |>
      filter(p.adjusted < 0.05) |>
      pluck("ID")
  })


go_data <- read_tsv(M$go_reference)
ontologizer <- get_ontologizer(M$ontologizer_path)
all_ontologizer <- bind_rows(
  mutate(ontologizer$id_with_open, subset = "id with open"),
  mutate(ontologizer$unknown_to_db, subset = "not DBP")
) |>
  bind_rows(
    bind_rows(intensity_tbs),
    bind_rows(enriched_mods)
  ) |>
  inner_join(go_data, by = join_by(x$ID == y$GO_IDs)) |>
  filter(p.adjusted < 0.05)

intensity_vecs_ontology <- lapply(intensity_vecs, \(x) {
  ids_into_ontology(x, target = "GOID", collapse = FALSE)
})
TABLES$all_ontologizer_sig <- all_ontologizer

GRAPHS$intensity_overlap <- ggVennDiagram(intensity_vecs) + scale_fill_paletteer_c("ggthemes::Classic Red")

# ----------------------------------------
# Verification for new grouping strategy

# Show that lfq intensity is consistent when grouping by unmatched peptides
# stderrs <- list()
# stds <- list()
# grouping_cols <- c("GroupUP", "Group")
# data <- inner_join(M$data, M$lfq, by = join_by(ProteinId))
# for (g in grouping_cols) {
#   by_intensity <- M$data |>
#     inner_join(merge_lfq(data, "mean"), by = join_by(ProteinId)) |>
#     select(ProteinId, log_intensity, {{ g }})
#   summarized <- by_intensity |>
#     filter(!is.na(log_intensity)) |>
#     group_by(!!as.symbol(g)) |>
#     summarize(intensity_std = sd(log_intensity, na.rm = TRUE))
#   stderrs[[g]] <- sd(summarized$intensity_std, na.rm = TRUE)
#   stds[[g]] <- summarized$intensity_std
# }

# stderr was lower when grouping by unmatched peptides, indicating this is a better way
# to form protein groups between different engines

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
