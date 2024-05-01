library(tidyverse)
library(glue)
args <- list(python_source = "./bin")

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

fixSep <- function(tb) {
  fixed <- tb %>% mutate(GO = map_chr(GOs,
                                      \(x) str_replace_all(x, ",", ";")))
}

filesToFix <- function(pattern) {
  if (str_detect(getwd(), "Bio_SDD")) {
    wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  } else {
    wd <- "/home/shannc/workflow"
  }
  return(list.files(glue("{wd}/results"), pattern = pattern, full.names = TRUE,
                    recursive = TRUE))
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
      tb %>% select(-Group) %>% write_tsv(., new_name)
     } else if (str_detect(filename, "intersected")) {
       tb %>% mutate(ID_method = "standard") %>% write_tsv(., filename)
     } else {
       tb %>% select(-Group) %>% write_tsv(., filename)
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
}

applyFixes <- function(file_list, fix_name) {
  lapply(file_list, \(x) fix(x, fix_name))
}


