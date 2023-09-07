library(tidyverse)
library(optparse)
library(Biostrings)



exact_coverage <- function(protein, peps) {
  # Based on pyteomics' implementation of coverage
  # is an exact calculation, no mismatches allowed
  protein <- gsub("[^A-Z]+", "", protein)
  mask <- rep(0, nchar(protein))
  all_locs <- lapply(peps, function(x) {
    return(str_locate_all(protein, x)[[1]])
  })
  for (l in all_locs) {
    if (dim(l)[1] == 0) {
      return(0)
    }
    sapply(seq_len(dim(l)[1]), function(x) {
      start <- l[x, "start"] %>% unlist(use.names = FALSE)
      mask[unname(l[x, "start"]):unname(l[x, "end"])] <<- 1
    })
  }
  return(sum(mask) / length(mask))
}

remove_substr <- function(pair) {
  # If any substring of peptide A is a substring of peptide B,
  # return peptide B and vice versa
  if (str_detect(pair[2], pair[1]) || str_detect(pair[1], pair[2])) {
    if (nchar(pair[2]) > nchar(pair[1])) {
      return(pair[2])
    } else {
      return(pair[1])
    }
  } else {
    return(pair)
  }
}

coverage <- function(protein, peps) {
  # Using alignment percent identity for the coverage
  # calculation takes into account gaps and mismatches
  if (any(protein == peps)) {
    return(1.00)
  }
  if (length(peps) > 1) {
    combos <- combn(peps, 2)
    peps <- unlist(lapply(seq_len(dim(combos)[2]), function(x) {
      return(remove_substr(combos[, x]))
    })) %>% unique()
  }
  align <- pairwiseAlignment(
    pattern = peps, subject = protein,
    type = "global-local"
  )
  cov <- Biostrings::coverage(align)
  cov <- data.frame(values = cov@values, lengths = cov@lengths)
  sum_cov <- cov %>% filter(values > 0) %>% select(lengths) %>% sum()
  return(sum_cov/nchar(protein))
}

clean_split_peps <- function(pep_list) {
  splits <- unlist(str_split(pep_list, ","), use.names = FALSE)
  splits <- splits[!(grepl("^NA$", splits))] %>% unique()
  cleaned <- unlist(lapply(splits, function(x) {
      return(paste0(str_extract_all(x, "[A-Z]+")[[1]], collapse = ""))
    }))
  return(cleaned)
}


coverage_calc <- function(prot_df) {
  covered_df <- prot_df %>%
    mutate(coverage = unlist(lapply(seq_along(prot_df$peptideIds), function(x) {
      peps <- clean_split_peps(.[x, ]$peptideIds)
      if (length(peps) > 1) {
        pairs <- combn(peps, 2)
        peps <- unlist(lapply(seq_len(dim(pairs)[2]), function(x) {
          remove_substr(pairs[, x])
        }), use.names = FALSE)
      }
      prot <- .[x, ]$sequence
      return(coverage(prot, peps))
    })))
  return(covered_df)
}


if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("--intersected_searches"),
    type = "character",
    help = "intersected_searches"
  )
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  args <- parse_args(parser)
  prot_df <- read.delim(args$intersected_searches, sep = "\t") %>% as_tibble()
  covered <- coverage_calc(prot_df)
  write_delim(covered, args$output, delim = "\t")
}

