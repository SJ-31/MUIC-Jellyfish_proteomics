entry_name_from_header <- function(header) {
  uniprot_regex <- "\\S*\\|\\S*\\|\\S* (.*) OS="
  ncbi_regex <- "\\S* (.*) \\["
  if (str_detect(header, ".*\\|.*\\|.*OS=.*")) {
    return(str_extract(header, uniprot_regex, group = 1))
  } else if (str_detect(header, "\\S* \\[.*\\].*")) {
    return(str_extract(header, ncbi_regex, group = 1))
  }
  header
}

ANNO_COLS <- c(
  "PFAMs", "header", "PANTHER", "eggNOG_description",
  "interpro_description", "interpro_accession", "eggNOG_OGs", "GO_IDs"
)
pfam_map <- read_tsv(glue("{M$wd}/data/reference/Pfam-A.clans.tsv"))
get_pfam_acc <- function(pfam_vec, pfam_map_file) {
  if (!exists("PFAM_SHORT_NAME2ACC")) {
    pfam_map <- read_tsv(pfam_map_file)
    short_name2acc <- setNames(as.list(pfam_map$accession), pfam_map$short_name)
    assign("PFAM_SHORT_NAME2ACC", short_name2acc, envir = globalenv())
  }
  map_chr(pfam_vec, \(x) {
    if (is.na(x)) {
      x
    } else if (str_sub(x, 1, 2) == "PF" && str_detect(x, "_")) {
      str_extract(x, "(PF.*?)_", group = 1)
    } else if (x %in% names(short_name2acc)) {
      short_name2acc[[x]]
    } else {
      NA
    }
  })
}

str_split_1_ignore_na <- function(x, pattern) {
  if (is.na(x)) {
    x
  } else {
    str_split_1(x, pattern)
  }
}

#' Various transformations to prepare `tb` to have
#' its entries categorized
#'
#' @description
#'
format_to_categorize <- function(tb, pfam_map_file, grouping_col) {
  tb <- select(tb, all_of("ProteinId", c(ANNO_COLS))) |> mutate(
    name = map_chr(header, entry_name_from_header),
    PFAM_accession = lapply(PFAMs, \(x) {
      str_split_1_ignore_na(x, ";") |>
        get_pfam_acc(pfam_map_file = pfam_map_file) |>
        unique()
    }),
    # TODO: check the others
  )
  tb |>
    group_by(!!as.symbol(grouping_col)) |>
    summarise(across(everything(), list)) # Group into list to make
  # element access easier (don't need to rely on string matching)
}
