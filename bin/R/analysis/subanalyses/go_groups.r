if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
GRAPHS <- list()
TABLES <- list()

go_data <- read_tsv(M$go_reference)
grouped_cat <- read_tsv(glue("{M$chosen_path}/Analysis/grouped_annotations_only.tsv"))

with_go <- function(tb) {
  separate_longer_delim(tb, GO_IDs, ";") |>
    group_by(GO_IDs) |>
    summarise(count = n(), assigned_COG = dplyr::first(assigned_COG)) |>
    inner_join(go_data, by = join_by(GO_IDs)) |>
    arrange(desc(count))
}

get_go <- function(category) {
  grouped_cat |>
    filter(assigned_COG == category) |>
    select(GO_IDs, assigned_COG) |>
    with_go()
}

GET_FILES <- FALSE
if (GET_FILES) {
  lapply(unique(grouped_cat$assigned_COG), \(x) {
    name <- str_replace_all(x, "/", "_") |> str_replace_all(" ", "_")
    TABLES[[glue("{name}")]] <<- get_go(x)
  })

  cog_names <- unique(grouped_cat$assigned_COG)
  cleaned_names <- cog_names |> map_chr(\(x) str_replace_all(x, "/", "_") |> str_replace_all(" ", "_"))

  key <- setNames(cog_names, cleaned_names)

  enriched_cog <- lapply(cleaned_names, \(x) {
    if (!is.na(x)) {
      tb <- read_tsv(glue("{M$ontologizer_path}/COG/{x}.tsv")) |> mutate(cog = key[[x]])
      TABLES[[glue("{x}_SIG_ENRICHED")]] <<- tb |>
        filter(p.adjusted < 0.05)
      tb
    }
  }) |>
    bind_rows()
}

# Singleton de novo peptide groups
cog_go_view <- function(tb, cog, go = TRUE) {
  if (go == TRUE) {
    filter(tb, assigned_COG == cog) |>
      select(GO_IDs, assigned_COG) |>
      with_go() |>
      see()
  } else {
    filter(tb, assigned_COG == cog) |> see()
  }
}

d <- grouped_cat |>
  filter(str_detect(header, "DENOVO"))
# cog_go_view(d, "Cytoskeleton", FALSE)
# cog_go_view(d, "Cytoskeleton")

fragment_names <- c("fragment", "partial")
fragment_regex <- paste0(fragment_names, collapse = "|")
frag <- grouped_cat |> filter(str_detect(str_to_lower(entry_name), fragment_regex))
# see(frag)

# save(c(GRAPHS, TABLES), glue("{M$outdir}/go_terms"))

joined_tox <- read_tsv(glue("{M$chosen_path}/Analysis/grouped_tox.tsv"))

with_go(joined_tox) |> see()
