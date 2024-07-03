if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

GRAPHS <- list()
TABLES <- list()

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


get_engine_counts <- function(type, run_data) {
  if (type == "open") {
    engines <- c("msfraggerGlyco", "msfraggerGPTMD", "metamorpheusGTPMD")
  } else {
    engines <- c("tide", "msfragger", "metamorpheus", "identipy", "comet", "msgf")
  }

  get_helper <- function(pass, run_tb) {
    if (type == "open") {
      dir <- glue("{M$path}/{pass}/Open_search/Percolator")
    } else {
      dir <- glue("{M$path}/{pass}/Percolator")
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

  passes <- list("1-First_pass", "2-Second_pass")
  metrics <- list()
  lapply(c(1, 2), \(x) {
    pass <- passes[[x]]
    metrics[[pass]] <- get_helper(pass, get_all_ids(run_data[[x]]))
  }) |> bind_rows()
}

get_mapped_scans <- function(pass) {
  list.files(glue("{M$path}/{pass}/Quantify/Mapped_scans"), pattern = "*.tsv", full.names = TRUE) |>
    lapply(read_tsv) |>
    bind_rows()
}

count_spectra <- function(
    spectra_filename, mapped_scans,
    spectra_path = glue("{M$wd}/data/MS/mzML")) {
  ms <- readMSData(glue("{spectra_path}/{spectra_filename}.mzML"), mode = "onDisk")
  header <- header(ms) |> as_tibble()
  spectra_levels <- header$msLevel |> table()
  current_scans <- mapped_scans |> filter(file == filename)
  ms1_current <- filter(current_scans, msLevel == 1)
  ms2_current <- filter(current_scans, msLevel == 2)
  ms_metrics <- tibble(
    n_ms1_spectra = spectra_levels[1],
    n_ms2_spectra = spectra_levels[2],
    n_psms_ms1 = nrow(ms1_current), # Number f
    n_matched_ms1 = length(unique(ms1_current$scan)),
    n_matched_ms2 = nrow(ms2_current), # Number of psms generated from ms2 spectra
    # a given spectra ms2 spectra can form multiple psms (saw at most 2)
    n_matched_ms2 = length(unique(ms2_current$scan)), # Number of ms2 spectra
    # that were matched, i.e. the unique ms2 spectra that matched
    n_unmatched =
    )
  engine_metrics <- current_scans$engine |>
    table() |>
    table2tb(id_col = "engine") |>
    mutate(file = filename) |>
    rename(n_matched_spectra = n)
  list(ms_metrics, engine_metrics)
}

get_unmatched_mzml <- function(pass) {
  list.files(glue("{M$path}/{pass}/Quantify/Unmatched"), pattern = "*.mzML") |>
    map_chr(\(x) str_replace(x, ".mzML", ""))
}


TABLES$standard_search_metrics <- get_engine_counts("standard", M$run)
TABLES$open_search_metrics <- get_engine_counts("open", M$run)

unmatched_ms <- readMSData(glue("{M$path}/{pass}/Quantify/Unmatched/filtered_{spectra_filename}.mzML"),
  mode = "onDisk"
)

mzmls <- list.files(glue("{M$wd}/data/MS/mzML"), pattern = "*.mzML") |> map_chr(\(x) str_replace(x, ".mzML", ""))

spectra_metrics <- lapply(mzmls, \(x) {
  count_spectra(x, get_mapped_scans("1-First_pass"))
}) |> purrr::reduce(
  \(x, y) list(bind_rows(x[[1]], y[[1]]), bind_rows(y[[1]], y[[2]]))
)



unmatched_mzmls_first <- get_unmatched_mzml("1-First_pass")



# save(TABLES, glue("{M$outdir}/search_metrics"))
