library(tidyverse)
library(glue)
args <- list(
  python_source = "./bin",
  go_path = "./data/reference/go.obo",
  go_slim_path = "./data/reference/goslim_generic.obo"
)


#' Add `engine` column in the unmatched peptide file showing which engines mapped
#' to that peptide
#'
find_engine_unmatched <- function(percolator_dir, unmatched_path) {
  get_percolator_peptides <- function() {
    transform_seq <- function(tb, engine) {
      tb |> mutate(
        peptide = map_chr(peptide, \(x) {
          str_replace_all(x, "[\\.-]", "")
        }),
        engine = engine
      )
    }
    file_list <- list.files(percolator_dir, pattern = "*_percolator_psms.tsv", full.names = TRUE)
    engines <- file_list %>%
      map_chr(., \(x) {
        gsub("_percolator_psms.tsv", "", x) %>%
          gsub(".*/", "", .)
      })
    tsvs <- lapply(file_list, read_tsv) %>% `names<-`(engines)
    if ("tide" %in% names(tsvs)) {
      tsvs$tide <- tsvs$tide |> rename(peptide = sequence)
    }
    lmap(tsvs, \(x) transform_seq(x[[1]], names(x))) %>% `names<-`(engines)
  }

  unmatched <- read_tsv(unmatched_path)
  peps <- get_percolator_peptides()
  to_join <- list()
  for (e in names(peps)) {
    to_join[[e]] <- inner_join(unmatched, peps[[e]], by = join_by(x$peptideIds == y$peptide)) |>
      select(ProteinId, engine)
  }
  all <- left_join(unmatched, bind_rows(to_join)) |>
    group_by(ProteinId) |>
    summarise(engine = paste0(engine, collapse = ";"))
  inner_join(unmatched, all, by = join_by(ProteinId))
}


#' Correct the "NCBI_ID" entry for database proteins that were originally derived
#' from UniProt
correctIds <- function(tb) {
  from_uniprot <- local({
    mapped <- map2_lgl(
      # Due to how the NCBI_ID column was generated, proteins from
      # UniProt would have the same entry in these two columns
      tb$UniProtKB_ID, tb$NCBI_ID,
      \(x, y) str_detect(x, y)
    )
    replace_na(mapped, FALSE)
  })
  tb[from_uniprot, ]$NCBI_ID <- NA
  return(tb)
}


# Temporary fixes for bugs in `sort_interpro` script
# - Removes duplicate ids from `matchedPeptideIds` column
# - Removes `query` column
cleanDuplicateIds <- function(tb) {
  if (!purrr::pluck_exists(tb, "ProteinId")) {
    ids <- tb$query
  } else {
    ids <- tb$ProteinId
  }
  cleaned_ids <- purrr::map2_chr(
    ids, tb$matchedPeptideIds,
    \(x, y) {
      if (!is.na(y) && grepl(pattern = x, y)) {
        cleaned <- str_split_1(y, ";") %>%
          discard(., \(z) z == x) %>%
          paste0(., collapse = ";")
        return(cleaned)
      } else {
        return(y)
      }
    }
  )
  tb$matchedPeptideIds <- cleaned_ids
  return(tb)
}

get_to_fix <- function(pattern) {
  if (str_detect(getwd(), "Bio_SDD")) {
    wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  } else {
    wd <- "/home/shannc/workflow"
  }
  return(list.files(glue("{wd}/results"),
    pattern = pattern, full.names = TRUE,
    recursive = TRUE
  ))
}

fix <- function(filename, fix) {
  tb <- read_tsv(filename)
  if (fix == "unify_groups") {
    # Sat Apr 20 18:01:36 2024
    # Fix for the error in `unify_groups`
    # Remove the "Group" column in "grouped_open_searches" files
    # Add "ID_method" in "intersected_searches" files
    if (grepl("grouped", filename)) {
      new_name <- gsub("grouped_", "", filename)
      tb %>%
        select(-Group) %>%
        write_tsv(., new_name)
    } else if (str_detect(filename, "intersected")) {
      tb %>%
        mutate(ID_method = "standard") %>%
        write_tsv(., filename)
    } else {
      tb %>%
        select(-Group) %>%
        write_tsv(., filename)
    }
  }
  if (fix == "sort_groups") {
    # Fri Apr 26 14:16:53 2024
    # Apply the correct grouping using the new unification method
    library(reticulate)
    source("./bin/R/combine_all.r")
    tb <- select(tb, -Group)
    tb <- unifyGroups(tb) %>% relocate(Group, .after = category)
    tb %>% write_tsv(., file = filename)
  }
  if (fix == "correct_ids") {
    # 2024-05-04 Fixed the identical entries in NCBI_ID and UniProtKB_ID
    source("./bin/R/combine_all.r")
    tb <- correctIds(tb)
    tb %>% write_tsv(., file = filename)
  }
  if (fix == "get_slims") {
    source("./bin/R/GO_helpers.r")
    tb <- tb %>%
      mutate(GO_slims = map_chr(GO_IDs, slimsFromGoString)) %>%
      relocate(., GO_slims, .after = GO_IDs)
    tb %>% write_tsv(., file = filename)
  }
  if (fix == "header_fix") {
    library("glue")
    # Wednesday 2024-06-06
    # Correct regex used to extract organism from header
    headerFix <- function(header) {
      return(str_extract(header, ".*\\[([A-Z].*)\\]", group = 1))
    }
    wrong <- tb |>
      filter(grepl("\\[", organism)) |>
      mutate(organism = map_chr(header, headerFix))
    print(glue("Dim wrong {dim(wrong)[1]}"))
    others <- tb %>% filter(!grepl("\\[", organism))
    print(glue("Dim others {dim(others)[1]}"))
    bind_rows(others, wrong) |> write_tsv(file = filename)
  }
  if (fix == "organism_header") {
    library("glue")
    source("./bin/R/helpers.r")
    # Friday 2024-06-14
    # Properly get organisms from header
    tb <- get_organism(tb)
    write_tsv(tb, file = filename)
  }
  if (fix == "group_unique_peps") {
    # 2024-06-23 Added a new, simpler method of grouping proteins
    # Just by their unique peptides
    source("./bin/R/helpers.r")
    groupByUniquePeptides(tb) |> write_tsv(file = filename)
  }
  if (fix == "mismatch") {
    rename(tb, n_mismatches = n_replacements) |> write_tsv(file = filename)
  }
}


apply_fixes <- function(file_list, fix_name) {
  lapply(file_list, \(x) fix(x, fix_name))
}
