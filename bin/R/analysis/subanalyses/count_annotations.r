if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

TABLES <- list()
# ----------------------------------------
# PFAM
pfam_ref <- read_tsv(glue("{M$path}/Databases/pfam_entries.tsv"))
pfams <- flatten_by(M$data$PFAMs, ";") |>
  table() |>
  table2tb(id_col = "from_sample")
has_ids <- pfams |> filter(map_lgl(
  pfams$from_sample,
  \(x) str_sub(x, 1, 2) == "PF" && str_detect(x, "_")
))

no_ids <- pfams |>
  filter(!from_sample %in% has_ids$from_sample) |>
  mutate(description = from_sample)

has_ids <- mutate(has_ids,
  id = map_chr(from_sample, \(x) str_extract(x, "(PF.*?)_", group = 1)),
  description = map_chr(id, \(x) {
    filtered <- pfam_ref |> filter(accession == x)
    if (nrow(filtered) == 0) {
      return(NA)
    }
    return(filtered$name)
  })
)

sample_pfams <- bind_rows(no_ids, has_ids) |> select(id, n, description)
TABLES$pfam <- sample_pfams |> mutate(database = "PFAM")


# ----------------------------------------
# PANTHER
panther_tb <- flatten_by(M$data$PANTHER, ";") |>
  table() |>
  table2tb(id_col = "from_sample") |>
  mutate(id = map_chr(from_sample, \(x) {
    str_extract(x, "(PTHR[0-9]+)[:_]", group = 1)
  }), description = map_chr(from_sample, \(x) {
    str_extract(x, "PTHR[0-9]+[:_](.*)", group = 1)
  }))

TABLES$panther <- panther_tb |>
  select(id, n, description) |>
  mutate(database = "PANTHER")

# ----------------------------------------
# Interpro
interpro_ref <- read_tsv(glue("{M$wd}/data/reference/interpro_entries.tsv"))
interpro_tb <- flatten_by(M$data$interpro_accession, ";") |>
  table() %>%
  table2tb(id_col = "id") |>
  left_join(interpro_ref, by = join_by(x$id == y$ENTRY_AC)) |>
  mutate(database = "interpro", description = paste0(ENTRY_TYPE, "| ", ENTRY_NAME)) |>
  select(id, n, description, database)

TABLES$interpro <- interpro_tb

# ----------------------------------------
# EGGNOG
eggnog_columns <- c("taxID", "OG", "COG_category", "description")
eggnog_ref <- read_tsv(glue("{M$wd}/data/reference/e5.og_annotations.tsv"),
  col_names = eggnog_columns
)

eggnog_tb <- flatten_by(M$data$eggNOG_OGs, ",") |>
  table() %>%
  table2tb(id_col = "id") |>
  mutate(id = map_chr(id, \(x) str_extract(x, "(.*)@", group = 1))) |>
  left_join(eggnog_ref, by = join_by(x$id == y$OG))

TABLES$eggnog <- eggnog_tb |>
  mutate(description = paste0("COG: ", COG_category, "|", description)) |>
  select(id, n, description) |>
  mutate(database = "eggNOG")

# ----------------------------------------
# GO
process_go <- function(tb, col) {
  proc <- flatten_by(tb[[col]], ";") |>
    table() |>
    table2tb(id_col = "id") |>
    left_join(go_ref, by = join_by(x$id == y$GO_IDs)) |>
    mutate(
      description = paste0(term, " (", ontology, ")", ": ", definition),
      database = "GO"
    ) |>
    select(id, n, description, database)
}

go_ref <- read_tsv(glue("{M$wd}/data/reference/go_data.tsv"))
go_tb <- process_go(M$data, "GO_IDs")
slim_tb <- process_go(M$data, "GO_slims")

TABLES$go_slims <- slim_tb
TABLES$go <- go_tb

save(c(TABLES), glue("{M$outdir}/Annotation_counts"))
