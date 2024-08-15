library("paletteer")
library("ggVennDiagram")
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
TABLES <- list()
GRAPHS <- list()
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
    purrr::pluck("GroupUP")
}



data <- M$data |>
  inner_join(M$lfq, by = join_by(ProteinId)) |>
  distinct(ProteinId, .keep_all = TRUE)
if (!file.exists(glue("{M$ontologizer_path}/high_intensity.tsv"))) {
  ont <- new.env()
  reticulate::source_python(glue("{M$python_source}/ontologizer_wrapper.py"), envir = ont)
  by_intensity <- merge_lfq(data, "mean") %>%
    inner_join(., dplyr::select(data, ProteinId, GO_IDs, GroupUP), by = join_by(ProteinId)) |>
    filter(!is.na(log_intensity)) |>
    group_by(GroupUP) |>
    summarise(GO_IDs = paste0(GO_IDs, collapse = ";"), log_intensity = mean(log_intensity)) |>
    mutate(GO_IDs = map_chr(GO_IDs, split_unique_join)) |>
    filter(!is.na(GO_IDs))
  O <- ont$Ontologizer(by_intensity, M$ontologizer_exec, M$go_path, "GroupUP")
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
  grouped <- data |>
    group_by(GroupUP) |>
    summarise(GO_IDs = paste0(GO_IDs, collapse = ";")) |>
    mutate(GO_IDs = map_chr(GO_IDs, split_unique_join)) |>
    filter(!is.na(GO_IDs))

  has_mods <- read_tsv(M$percolator_all_path) |>
    inner_join(data, by = join_by(ProteinId)) |>
    filter(GroupUP %in% grouped$GroupUP) |>
    filter(!is.na(mods)) |>
    fix_all_mods()

  all_mods <- flatten_by(has_mods$mods, ";")
  mod_table <- all_mods |>
    map_chr(\(x) str_extract(x, "(.*\\|.*)\\|[1-9]+", group = 1)) |>
    discard(\(x) str_detect(x, "229.16")) |> # This shouldn't be here...
    table()

  has_mods_ids <- lapply(c("Met\\|15.99", "Nterm\\|42.0", "Lys\\|42.01"), \(x) {
    filter(has_mods, grepl(x, mods)) |>
      pluck("GroupUP") |>
      unique()
  }) |>
    `names<-`(mod_names)
  ont <- new.env()

  reticulate::source_python(glue("{M$python_source}/ontologizer_wrapper.py"), envir = ont)
  O <- ont$Ontologizer(grouped, M$ontologizer_exec, M$go_path, id_col = "GroupUP")
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


all_ontologizer$subset |> table()
see(all_ontologizer)


slims <- get_go_slim(unique(all_ontologizer$ID), M$go_path, M$go_slim_path)
all_ontologizer$slim <- map_chr(all_ontologizer$ID, \(x) {
  find_slims <- slims[[x]]$all |> discard(\(y) y == x)
  if (length(find_slims) > 0) {
    paste0(find_slims, collapse = ";")
  } else {
    "NONE"
  }
})

all_ontologizer |> see()

TABLES$all_ontologizer_sig <- all_ontologizer



intensity_vecs_ontology
all_ontologizer$ontology

GRAPHS$intensity_overlap <- ggVennDiagram(intensity_vecs) + scale_fill_paletteer_c("ggthemes::Classic Red")

save(c(TABLES, GRAPHS), glue("{M$outdir}/misc"))
