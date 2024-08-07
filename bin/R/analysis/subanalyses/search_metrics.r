options(
  browser = "firefox",
  rlang_backtrace_on_error = "full",
  error = rlang::entrace
)
rlang::global_entrace()
library("MSnbase")

if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

GRAPHS <- list()
TABLES <- list()
PASSES <- list("1-First_pass", "2-Second_pass")

get_engine_stats <- function(engine, percolator_dir, all_protein_ids, fdr = 0.05) {
  protein_file <- glue("{percolator_dir}/{engine}_percolator_proteins.tsv")
  if (!file.exists(protein_file)) {
    return(tibble())
  }
  protein_groups <- read_tsv(protein_file)
  proteins <- protein_groups |> separate_longer_delim("ProteinId", ",")
  psms <- read_tsv(glue("{percolator_dir}/{engine}_percolator_psms.tsv"))
  peptides <- read_tsv(glue("{percolator_dir}/{engine}_percolator_peptides.tsv"))
  if (engine == "tide") {
    psms <- psms |> rename(`q-value` = `percolator q-value`, proteinIds = `protein id`)
    peptides <- peptides |> rename(`q-value` = `percolator q-value`, proteinIds = `protein id`)
  }
  row <- tibble(
    engine = engine,
    n_groups = nrow(protein_groups),
    n_groups_with_duplicates = protein_groups |> filter(grepl(",", ProteinId)) |> nrow(),
    n_proteins = nrow(proteins),
    n_groups_below_fdr = protein_groups |> filter(`q-value` <= fdr) |> nrow(),
    n_proteins_below_fdr = proteins |> filter(`q-value` <= fdr) |> nrow(),
    n_proteins_accepted = proteins |> filter(ProteinId %in% all_protein_ids) |> nrow(),
    # To be accepted, a protein needs to be matched by at least two other engines
    # And the q value of at least two psm matches needs to be below the fdr
    n_denovo_proteins = proteins |> filter(grepl("D", ProteinId)) |> nrow(),
    n_database_proteins = proteins |> filter(grepl("P", ProteinId)) |> nrow(),
    n_transcriptome_proteins = proteins |> filter(grepl("T", ProteinId)) |> nrow(),
    n_psms = nrow(psms),
    n_psms_below_fdr = psms |> filter(`q-value` <= fdr) |> nrow(),
    n_unmatched_psms = filter(psms, is.na(nchar(proteinIds))) |> nrow(),
    n_peptides = nrow(peptides),
    n_peptides_below_fdr = peptides |> filter(`q-value` <= fdr) |> nrow(),
    n_unmatched_peptides = filter(peptides, is.na(nchar(proteinIds))) |> nrow()
  )
  row
}

get_all_ids <- function(tb) {
  c(
    tb$ProteinId,
    unlist(lapply(discard(tb$MatchedPeptideIds, is.na), str_split_1, pattern = ";"))
  )
}


get_engine_counts <- function(type, run_data, path = M$path) {
  if (type == "open") {
    engines <- c("msfraggerGlyco", "msfraggerGPTMD", "metamorpheusGTPMD")
  } else {
    engines <- c("tide", "msfragger", "metamorpheus", "identipy", "comet", "msgf")
  }

  get_helper <- function(pass, run_tb) {
    if (type == "open") {
      dir <- glue("{path}/{pass}/Open_search/Percolator")
    } else {
      dir <- glue("{path}/{pass}/Percolator")
    }
    metrics <- lapply(engines, \(x) {
      get_engine_stats(x,
        percolator_dir = dir,
        all_protein_ids = run_tb
      )
    }) |>
      bind_rows() |>
      mutate(pass = pass)
  }

  metrics <- list()
  lapply(c(1, 2), \(x) {
    pass <- PASSES[[x]]
    metrics[[pass]] <- get_helper(pass, get_all_ids(run_data[[x]]))
  }) |> bind_rows()
}

get_mapped_scans <- function(pass, mode = "standard", path = M$path) {
  if (mode == "standard") {
    path <- glue("{path}/{pass}/Quantify/Mapped_scans")
  } else {
    path <- glue("{path}/{pass}/Open_search/Mapped_scans")
  }
  list.files(path, pattern = "*.tsv", full.names = TRUE) |>
    lapply(read_tsv) |>
    bind_rows()
}

plot_param_psms <- function(file_list, palette) {
  tb <- lapply(file_list, \(x) {
    param <- param_name_from_file(x)
    read_tsv(x) |>
      mutate(param = param) |>
      mutate(file = map_chr(file, \(x) str_remove_all(x, "msc_|filtered_|-calib")))
  }) |>
    bind_rows()

  longer <- tb |>
    pivot_longer(cols = where(is.numeric)) |>
    group_by(name, param, pass) |>
    summarize(value = sum(value)) |>
    ungroup() |>
    rename(Engine = name, `PSM count` = value, Pass = pass) |>
    rename_passes()

  plot_grouped_bar(longer, "PSM count", palette = palette)
}

count_spectra <- function(spectra_path, mapped_scans) {
  filename <- str_extract(spectra_path, ".*/(.*).mzML", group = 1)
  ms <- readMSData(spectra_path, mode = "onDisk")
  header <- header(ms) |> as_tibble()
  spectra_levels <- header$msLevel |> table()
  current_scans <- mapped_scans |> filter(file == filename)
  ms1_current <- filter(current_scans, msLevel == 1)
  ms2_current <- filter(current_scans, msLevel == 2)
  ms_metrics <- tibble(
    file = filename,
    n_ms1_spectra = spectra_levels[1],
    n_ms2_spectra = spectra_levels[2],
    n_psms_ms1 = nrow(ms1_current),
    n_matched_ms1 = length(unique(ms1_current$scan)),
    n_psms_ms2 = nrow(ms2_current), # Number of psms generated from ms2 spectra
    # a given spectra ms2 spectra can form multiple psms (saw at most 2)
    n_matched_ms2 = length(unique(ms2_current$scan)), # Number of ms2 spectra
    # that were matched, i.e. the unique ms2 spectra that matched
  )
  engine_metrics <- current_scans$engine |>
    table() |>
    table2tb(id_col = "engine") |>
    mutate(file = filename) |>
    rename(n_matched_spectra = n)
  list(ms_metrics, engine_metrics)
}

get_unmatched_mzml <- function(pass, path = M$path) {
  list.files(glue("{path}/{pass}/Quantify/Unmatched"), pattern = "*.mzML", full.names = TRUE)
}

main_metrics <- function(spectra_file_list, pass, unmatched_ms = FALSE, path = M$path) {
  if (unmatched_ms) {
    scans <- get_mapped_scans(pass, mode = "open", path = path)
  } else {
    scans <- get_mapped_scans(pass, path = path)
  }
  lst <- lapply(spectra_file_list, \(x) {
    count_spectra(x, scans)
  })
  if (length(lst) == 0) {
    stop("lst empty")
  }
  reduced <- lst |> purrr::reduce(
    \(x, y) list(bind_rows(x[[1]], y[[1]]), bind_rows(x[[2]], y[[2]]))
  )
  engines <- reduced[[2]] |>
    pivot_wider(names_from = engine, values_from = n_matched_spectra) |>
    mutate(pass = pass) # Number of engine psms for both ms1 and ms2 spectra
  spectra <- reduced[[1]] |> mutate(pass = pass)
  return(list(engine_n_psms_12 = engines, spectra_metrics = spectra))
}

param_name_from_file <- function(filename) {
  str_remove(filename, ".*/") |> str_remove("_.*.tsv")
}

rename_passes <- function(tb) {
  pass_lab <- get_label_replacement(c("First" = "1-First_pass", "Second" = "2-Second_pass"))
  tb |>
    mutate(Pass = pass_lab(Pass))
}

plot_grouped_bar <- function(tb, show_col, palette) {
  ggplot(tb, aes(x = Engine, pattern = Pass, y = !!as.symbol(show_col), fill = Engine)) +
    geom_bar_pattern(
      position = "dodge", stat = "identity",
      color = "black", pattern_density = 0.1, pattern_spacing = 0.03
    ) +
    scale_pattern_manual(values = c(First = "none", Second = "stripe")) +
    theme(
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 15), legend.text = element_text(size = 13),
      strip.text.x = element_text(size = 16), legend.title = element_text(size = 14)
    ) +
    facet_wrap(~param) +
    scale_fill_paletteer_d(palette)
}

mzmls <- list.files(glue("{M$wd}/data/MS/mzML"), pattern = "*.mzML", full.names = TRUE)
calib_mzmls <- list.files(glue("{M$path}/Preprocessed/Metamorpheus"), pattern = "*.mzML", full.names = TRUE)
pw_mzmls <- list.files(glue("{M$path}/Preprocessed/Proteowizard"), pattern = "*.mzML", full.names = TRUE)


run_paths <- c(M$path, M$cpath, M$mpath, M$ndpath)
run_names <- c("default", "Calibrated", "msConvert", "ND")
prefixes <- c("C_indra", "C_indra.calibrated", "C_indra.msconvert", "ND_C_indra")
all_mzmls <- list(mzmls, calib_mzmls, pw_mzmls, mzmls)

previous_saved <- list.files(glue("{M$wd}/docs/figures/search_metrics"), pattern = ".*tsv", full.names = TRUE)

if (all(map_lgl(run_names, \(x) any(grepl(x, previous_saved))))) {
  GET_SPECTRA <- TRUE
} else {
  GET_SPECTRA <- FALSE
}

if (!GET_SPECTRA) {
  for (i in seq_along(all_mzmls)) {
    metrics_1 <- main_metrics(
      all_mzmls[[i]],
      PASSES[[1]],
      path = run_paths[i]
    )
    metrics_2 <- main_metrics(
      all_mzmls[[i]],
      PASSES[[2]],
      path = run_paths[i]
    )
    unmatched_metrics_1 <- main_metrics(
      get_unmatched_mzml(PASSES[[1]], run_paths[i]), PASSES[[1]], TRUE, run_paths[i]
    )
    unmatched_metrics_2 <- main_metrics(
      get_unmatched_mzml(PASSES[[2]], run_paths[i]), PASSES[[2]], TRUE, run_paths[i]
    )

    TABLES[[glue("{run_names[i]}_spectra_metrics")]] <- bind_rows(
      metrics_1$spectra_metrics,
      metrics_2$spectra_metrics
    )
    TABLES[[glue("{run_names[i]}_engine_psm_counts")]] <- bind_rows(
      metrics_1$engine_n_psms_12,
      metrics_2$engine_n_psms_12
    )
    TABLES[[glue("{run_names[i]}_open_engine_psm_counts")]] <- bind_rows(
      unmatched_metrics_1$engine_n_psms_12,
      unmatched_metrics_2$engine_n_psms_12
    )
    TABLES[[glue("{run_names[i]}_unmatched_spectra_metrics")]] <- bind_rows(
      unmatched_metrics_1$spectra_metrics,
      unmatched_metrics_2$spectra_metrics
    )
    cur_run <- get_run(prefixes[i], run_paths[i])
    TABLES[[glue("{run_names[i]}_standard_search_metrics")]] <- get_engine_counts("standard", cur_run, run_paths[i])
    TABLES[[glue("{run_names[i]}_open_search_metrics")]] <- get_engine_counts("open", cur_run, run_paths[i])
  }
  save(c(TABLES, GRAPHS), glue("{M$wd}/docs/figures/search_metrics"))
  q()
}
library("ggpattern")

standard_search <- keep(previous_saved, \(x) str_detect(x, "engine_psm_counts") & !str_detect(x, "open"))
open_search <- keep(previous_saved, \(x) str_detect(x, "engine_psm_counts") & str_detect(x, "open"))
GRAPHS$standard_psms <- plot_param_psms(standard_search, "ggthemes::excel_Depth")
attr(GRAPHS$standard_psms, "height") <- 5
GRAPHS$open_psms <- plot_param_psms(open_search, "ggthemes::Classic_Purple_Gray_12")
attr(GRAPHS$open_psms, "height") <- 5


palettes <- c("ggthemes::excel_Vapor_Trail", "ggthemes::hc_darkunica")
os <- c("open", "standard")
for (i in seq_along(os)) {
  group_files <- keep(previous_saved, \(x) str_detect(x, glue("{os[i]}_search_metrics")))
  search_metrics <- lapply(group_files, \(x) {
    param <- param_name_from_file(x)
    read_tsv(x) |>
      mutate(param = param) |>
      rename(Engine = engine, Pass = pass) |>
      rename_passes()
  }) |>
    bind_rows() |>
    mutate(n_groups = log2(n_groups))


  show_col <- "n_groups"
  GRAPHS[[glue("n_groups_{os[i]}")]] <- plot_grouped_bar(
    search_metrics, show_col,
    palettes[i]
  ) + ylab("log2 n groups")
  if (os[i] == "open") {
    GRAPHS[[glue("n_groups_{os[i]}")]] <- GRAPHS[[glue("n_groups_{os[i]}")]] + theme(axis.title.y = element_blank())
  }
  groups_below <- lapply(run_names, \(x) {
    f <- filter(search_metrics, param == x)
    f[[show_col]]
  }) |> `names<-`(run_names)

  test_groups_below <- kruskal.test(groups_below)
  if (test_groups_below$p.value > 0.05) {
    TABLES[[glue("kruskal_groups_{os[i]}_not_sig")]] <- 0
  }
}

spectra_files <- keep(previous_saved, \(x) {
  str_detect(x, "spectra_metrics") & !str_detect(x, "unmatched") & !str_detect(x, "ND")
})

TABLES$standard_spectra <- spectra_files |>
  lapply(\(x) {
    p <- param_name_from_file(x)
    read_tsv(x) |>
      mutate(file = map_chr(file, \(x) str_remove_all(x, "msc_|filtered_|-calib"))) |>
      filter(pass == "1-First_pass") |>
      select(file, n_ms2_spectra) |>
      mutate(param = p)
  }) |>
  bind_rows() |>
  pivot_wider(id_cols = "file", names_from = "param", values_from = n_ms2_spectra) |>
  gt()


open_spectra_files <- keep(previous_saved, \(x) {
  str_detect(x, "spectra_metrics") & str_detect(x, "unmatched") & !str_detect(x, "ND")
})
open_spectra <- open_spectra_files |>
  lapply(\(x) {
    p <- param_name_from_file(x)
    read_tsv(x) |>
      mutate(file = map_chr(file, \(x) str_remove_all(x, "msc_|filtered_|-calib")), param = p) |>
      rename(Pass = pass) |>
      select(file, n_ms2_spectra, Pass, param) |>
      rename_passes()
  }) |>
  bind_rows()

GRAPHS$open_spectra_stats <- open_spectra |> ggplot(aes(x = file, y = n_ms2_spectra, fill = file, pattern = Pass)) +
  geom_bar_pattern(
    position = "dodge", stat = "identity",
    color = "black", pattern_density = 0.1, pattern_spacing = 0.03
  ) +
  scale_pattern_manual(values = c(First = "none", Second = "stripe")) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  facet_wrap(~param) +
  scale_fill_paletteer_d("lisa::BridgetRiley")
attr(GRAPHS$open_spectra_stats, "height") <- 6


save(c(TABLES, GRAPHS), glue("{M$wd}/docs/figures/search_metrics"))
