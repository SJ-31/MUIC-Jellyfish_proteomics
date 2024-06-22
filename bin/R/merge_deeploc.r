library("typed")
library("glue")
library("hash")
library("tidyverse")

getDeeploc <- function(deeploc_path, unmatched_path) {
  dl <- read_csv(deeploc_path)
  um <- read_tsv(unmatched_path)
  merged <- inner_join(um, dl,
    by = join_by(x$ProteinId == y$Protein_ID)
  ) %>%
    mutate(
      inferred_by = "deeploc",
      localization = Localizations,
      category = Signals,
    ) %>%
    select(all_of(c(colnames(um), "inferred_by", "localization", "category")))
  return(merged)
}


# Map localization data from deeploc onto GO
LOCALIZATION_MAP <- NULL
localization2go <- Character() ? function(loc = ? Character()) {
  if (is.null(LOCALIZATION_MAP)) {
    map <- hash::hash(
      "Cell membrane" = "GO:0005886",
      "Extracellular" = "GO:0005615",
      "Cytoplasm" = "GO:0005737",
      "Nucleus" = "GO:0005634",
      "Plastid" = "GO:0009536"
    )
    assign("LOCALIZATION_MAP", map, envir = rlang::global_env())
  }
  str_split_1(loc, "\\|") |>
    map_chr(\(x) LOCALIZATION_MAP[[x]]) |>
    paste0(collapse = ";")
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-m", "--main_results"), type = "character", action = "store")
  parser <- add_option(parser, c("-d", "--deeploc"), type = "character", action = "store")
  parser <- add_option(parser, c("-u", "--unmatched"), type = "character", action = "store")
  parser <- add_option(parser, c("-o", "--output"), type = "character", action = "store")
  args <- parse_args(parser)
  deeploc <- getDeeploc(args$deeploc, args$unmatched) |>
    mutate(GO_IDs = map_chr(localization, localization2go)) |>
    select(-localization)
  data <- read_tsv(args$main_results) |> bind_rows(deeploc)
  write_tsv(data, args$output)
}
