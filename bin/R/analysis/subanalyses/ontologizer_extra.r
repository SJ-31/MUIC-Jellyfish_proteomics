library("paletteer")
library("ggVennDiagram")
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
source(glue("{args$r_source}/GO_text_mining_helpers.r"))
source(glue("{args$r_source}/GO_chord.r"))
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
}

# ----------------------------------------
# Enrich terms based on modifications
mod_names <- c("Met_ox", "Nterm_acetyl", "Lys_acetyl")
if (!file.exists(glue("{M$ontologizer_path}/Met_ox.tsv"))) {
  grouped <- data |>
    group_by(GroupUP) |>
    summarise(GO_IDs = paste0(GO_IDs, collapse = ";")) |>
    mutate(GO_IDs = map_chr(GO_IDs, split_unique_join)) |>
    filter(!is.na(GO_IDs))

  hh <- new.env()
  reticulate::source_python(glue("{args$python_source}/helpers.py"), envir = hh)
  has_mods <- data |>
    filter(grepl("\\[", peptideIds)) |>
    separate_longer_delim(peptideIds, ";") |>
    mutate(mods = map_chr(peptideIds, hh$get_mods)) |>
    filter(!is.na(mods) & mods != "NA" & mods != "") |>
    fix_all_mods()


  write_tsv(has_mods, glue("{M$ontologizer_path}/has_mods.tsv"))


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

  mod_binary <- data |>
    select(GroupUP) |>
    distinct()
  for (s in names(has_mods_ids)) {
    mod_binary[[s]] <- mod_binary$GroupUP %in% has_mods_ids[[s]]
  }
  write_tsv(mod_binary, glue("{M$ontologizer_path}/mod_binary.tsv"))

  ont <- new.env()
  reticulate::source_python(glue("{M$python_source}/ontologizer_wrapper.py"), envir = ont)
  O <- ont$Ontologizer(grouped, M$ontologizer_exec, M$go_path, id_col = "GroupUP")
  params <- list(`-m` = "Bonferroni-Holm")

  enriched_mods <- O$runAll(has_mods_ids, params)
  lmap(enriched_mods, \(x) {
    write_tsv(x[[1]], glue("{M$ontologizer_path}/{names(x)}.tsv"))
  })
  enriched_mods <- lapply(enriched_mods, as_tibble)
}

if (!dir.exists(glue("{M$ontologizer_path}/COG"))) {
  dir.create(glue("{M$ontologizer_path}/COG"))
  ont <- new.env()
  reticulate::source_python(glue("{M$python_source}/ontologizer_wrapper.py"), envir = ont)
  data <- read_tsv(M$data_w_cat_path)
  with_cog <- group_by(data, GroupUP) |>
    summarise(
      GO_IDs = paste0(GO_IDs, collapse = ";")
    ) |>
    mutate(
      GO_IDs = map_chr(GO_IDs, split_unique_join),
    ) |>
    filter(!is.na(GO_IDs))
  O <- ont$Ontologizer(with_cog, M$ontologizer_exec, M$go_path, "GroupUP")
  cog_names <- unique(data$assigned_COG)
  cleaned_names <- cog_names |> map_chr(\(x) str_replace_all(x, "/", "_") |> str_replace_all(" ", "_"))
  groups <- lapply(cog_names, \(x) {
    data |>
      filter(assigned_COG == x) |>
      pluck("GroupUP") |>
      unique()
  }) |>
    `names<-`(cleaned_names)

  params <- list(`-m` = "Bonferroni-Holm")
  enriched_cog <- O$runAll(groups, params)

  lmap(enriched_cog, \(x) {
    write_tsv(x[[1]], glue("{M$ontologizer_path}/COG/{names(x)}.tsv"))
  })
}
