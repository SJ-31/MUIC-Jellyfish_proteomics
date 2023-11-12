library(hash)
library(tidyverse)
TEST <- TRUE

if (TEST) {
  test <- ("../results/jellyfish/1-First_pass/Combined/unified_groups.tsv")
  ms <- "../results/jellyfish/1-First_pass/Open_search/Percolator/msfraggerGPTMD_percolator_proteins.tsv" %>% read_tsv()
  mm <- "../results/jellyfish/1-First_pass/Open_search/Percolator/metamorpheusGTPMD_percolator_proteins.tsv" %>% read_tsv()
}

# CA = Common Artifact
# CB = Common Biological
# M = Metal
MOD_LIST_DEFAULT <- list("Cys_Carbamidomethylation" = "\\[57.02[0-9]*\\]",
                 "NTerm_Acetylation" = "\\[42.0[0-9]*\\]",
                 "Met_Oxidation" = "\\[15.99[0-9]*\\]")

cleanPeptide <- function(modified_pep) {
  clean_pep <- str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
    paste0(collapse = "")
  return(clean_pep)
}

namedModList <- function(df) {
  # Capture named modifications from dataframe as
  # a list to add to MOD_LIST
  capture_regex <- "^([A-Z])[a-z]+( ([A-Z])[a-z]+)?(:.*)"
  replacement_regex <- "\\1\\3\\4"
  mod_string <- "\\[([0-9]+\\.[0-9]+)|([a-zA-Z ]+:[ a-zA-Z]+)\\]"
  named_mod_string <- "[a-zA-Z ]+:[ a-zA-Z]+"
  unique_mods <- df$peptideIds %>%
    lapply(., str_extract, pattern = mod_string) %>%
    unlist(use.names = FALSE) %>%
    unique()
  named_mods <- str_extract(unique_mods, named_mod_string) %>%
    unique() %>%
    dropNA()
  if (length(named_mods) == 0) {
    return(list())
  }
  mod_names <- named_mods %>%
    str_replace(., capture_regex, replacement_regex) %>%
    str_replace_all(., " ", "_")
  named_mods <- paste0("\\[", named_mods, "\\]")
  named_mods <- as.list(named_mods) %>% `names<-`(mod_names)
  return(named_mods)
}

countMods <- function(peptide, mod_dict) {
  # Count the number of modifications in "peptide", then add to
  # the dictionary
  for (n in seq_along(MOD_LIST_DEFAULT)) {
    current_mod <- names(MOD_LIST_DEFAULT)[n]
    if (has.key(current_mod, mod_dict)) {
      current_val <- mod_dict[[current_mod]]
    } else {
      current_val <- 0
    }
    count <- grepl(MOD_LIST_DEFAULT[[n]], peptide)
    if (count > 0) {
        mod_dict[current_mod] <- count + current_val
    }
  }
}

getAllMods <- function(mod_df, row) {
  # Count all the modifications in a given peptide, returning
  # a string with the count of each mod separated by |
  current <- mod_df[row, ]
  mods <- hash()
  modified <- current$peptideIds %>%
  strsplit(",") %>%
  unlist() %>%
  unique() %>%
  grep("[][]", ., value = TRUE)
  lapply(modified, countMods, mod_dict = mods)
  mod_string <- paste(keys(mods), values(mods), collapse = "|")
  return(mod_string)
}


if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"), action = "store_true",
                       defaulta = TRUE, help = "Print extra output [default]")
  parser <- add_option(parser, c("-o", "--output"), type = "character",
                       help = "Output file name")
  args <- parse_args(parser)
  if (TEST) {
    args$input <- test
  }
  df <- read_tsv(args$input)
  has_mods <- df %>% filter(grepl("[][]", df$peptideIds))
  no_mods <- df[!(df$ProteinId %in% has_mods$ProteinId), ] %>%
    mutate(Mods = NA)
  MOD_LIST_DEFAULT <- c(MOD_LIST_DEFAULT, namedModList(has_mods))
  mod_col <- seq_len(dim(has_mods[1])) %>%
    lapply(., getAllMods, mod_df = has_mods) %>%
    unlist(use.names = FALSE)
  has_mods <- has_mods %>% mutate(Mods = mod_col)
  df <- bind_rows(has_mods, no_mods)
}
