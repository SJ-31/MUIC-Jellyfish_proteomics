
TEST <- FALSE

joinMods <- function(peptides) {
 pattern <- "([A-Z\\]]) ([A-Z\\[])"
 splits <- str_replace_all(peptides, pattern, "\\1,\\2") %>%
   unlist(use.names = FALSE)
 return(paste0(splits, collapse = ","))
}

dropNA <- function(vector) {
  return(vector[!is.na(vector)])
}

cleanNA <- function(vector) {
  cleaned <- lapply(vector, function(x) {
    if (is.na(x)) {
      return("")
    } else return(x)
  })
  return(unlist(cleaned, use.names = FALSE))
}

sortDuplicates <- function(df) {
  dupes <- df %>% filter(grepl(",", df$ProteinId))
  if (dim(dupes)[1] == 0) {
    return(df)
  }
  no_dupes <- df[!(df$ProteinId %in% dupes$ProteinId), ]
  split <- lapply(seq_len(dim(dupes)[1]), splitDuplicates,
                  dupe_table = dupes) %>% bind_rows()
  return(bind_rows(split, no_dupes))
}

cleanUp <- function(df, q_threshold) {
  # Format protein rows with multiple peptides (duplicates)
  # Separate duplicates into separate rows
  # Filter by q-value
  df <- mutate(df, peptideIds = unlist(lapply(df$peptideIds, joinMods),
                                       use.names = FALSE))
  df <- sortDuplicates(df) %>% filter(`q-value` <= 0.05)
  return(df)
}



if (TEST) {
  setwd('/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/')
  mm <- "../results/jellyfish/1-First_pass/Open_search/Percolator/metamorpheusGTPMD_percolator_proteins.tsv" %>% read_tsv()
  ms <- "../results/jellyfish/1-First_pass/Open_search/Percolator/msfraggerGPTMD_percolator_proteins.tsv" %>% read_tsv()
  ml <- "../results/jellyfish/1-First_pass/Open_search/Percolator/msfraggerGlyco_percolator_proteins.tsv" %>% read_tsv()
}


if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library(tidyverse)
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-v", "--verbose"), action = "store_true",
                       default = TRUE, help = "Print extra output [default]")
  parser <- add_option(parser, c("-i", "--input_directory"),
                       type = "character")
  parser <- add_option(parser, c("-r", "--r_source"),
                       type = "character")
  parser <- add_option(parser, c("-o", "--output"), type = "character",
                       help = "Output file name")
  parser <- add_option(parser, c("-m", "--seq_header_file"),
    type = "character",
    help = "Path to seq-header mapping"
  )
  args <- parse_args(parser)
  if (TEST) {
    args <- list(input_directory = "../results/jellyfish/1-First_pass/Open_search/Percolator/", seq_header_file = "../results/jellyfish/Databases/seq-header_mappings.tsv", r_source "~/Bio_SDD/MUIC_senior_project/workflow/bin")
  }
  setwd(args$r_source)
  source("./atleast2.r")
  mapped <- read_tsv(args$seq_header_file)
  tsv <- list.files(args$input_directory,
                    pattern = "*percolator_proteins.tsv") %>%
    paste0(args$input_directory, .) %>%
    lapply(., read_tsv)
  empty <- lapply(tsv, nrow)
  cleaned <- tsv[empty != 0] %>% lapply(., cleanUp)
  grouped <- bind_rows(cleaned) %>%
    mutate(Identification_method = "open_search") %>%
    group_by(ProteinId) %>%
    mutate_at(., vars(-group_cols()), paste0, collapse = ",") %>%
    distinct() %>%
    ungroup() %>%
    inner_join(mapped, by = join_by(x$ProteinId == y$id))
  write_delim(grouped, args$output, delim = "\t")
}
