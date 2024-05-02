library(tidyverse)
library(seqinr)
library(glue)
library(Biostrings)


writeAlignments <- function(row, file_name) {
  header <- ifelse(row[["header"]] == "unknown", row[["ProteinId"]],
    row[["header"]]
  )
  header <- glue("{row[['ProteinId']]}:{header}")
  pep_count <- row[["num_unique_peps"]]
  write.fasta(row[["seq"]], header, open = "a", file.out = file_name)
  write.fasta(row[["alignment"]],
    glue("ALIGNED PEPTIDES | COUNT: {pep_count}"),
    open = "a",
    file.out = file_name
  )
  cat("\n", file = file_name, append = TRUE)
}

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

splitWithGroups <- function(str) {
  group_count <- 0
  groups <- unlist(str_extract_all(str, "(\\[[A-Z]+\\])"))
  mask <- gsub("(\\[[A-Z]+\\])", "0", str)
  splits <- str_split_1(mask, "")
  replace_masks <- lapply(seq_along(splits), \(x) {
    if (splits[x] == "0") {
      group_count <- group_count + 1
      return(groups[group_count])
    }
    return(splits[x])
  }) %>% unlist()
  return(replace_masks)
}

#' Function for handling when identified peptides align to the same
#' indices on the query string
#'
#' @description
#' If the conflicting peptides have the same residues, then do nothing
#' if peptide1 has an "A" residue but peptide2 has a "C" residue, since
#' we don't know which is correct, this resolves into "[AC]"
resolveAlignment <- function(seq1, seq2) {
  if (grepl("\\[", seq1)) {
    split1 <- splitWithGroups(seq1)
  } else {
    split1 <- str_split_1(seq1, "")
  }
  split2 <- str_split_1(seq2, "")
  resolved <- sapply(seq_along(split1), \(x) {
    if (split1[x] == "-" && split2[x] == "-") {
      return("-")
    } else if (split1[x] != "-" && split2[x] == "-") {
      return(split1[x])
    } else if (split1[x] == "-" && split2[x] != "-") {
      return(split2[x])
    } else if (split1[x] == split2[x] || length(str_which(
      split1[x],
      split2[x]
    )) > 0
    ) {
      return(split1[x])
    } else if (grepl("\\[", split1[x])) {
      temp_split <- str_split_1(split1[x], "")
      return(paste0(R.utils::insert(temp_split, 2, split2[x]), collapse = ""))
    } else {
      return(paste0("[", split1[x], split2[x], "]", collapse = ""))
    }
  })
  return(str_flatten(resolved))
}


appendAllList <- function(lst, to_append) {
  if (length(to_append) > 1) {
    new <- lapply(seq_along(lst), \(x) {
      lst[[x]] <- c(lst[[x]], to_append[x])
    })
  } else {
    new <- lapply(seq_along(lst), \(x) {
      lst[[x]] <- c(lst[[x]], to_append)
    })
  }
  return(new)
}

coverage <- function(protein, peps) {
  # Using alignment percent identity for the coverage
  # calculation takes into account gaps and mismatches
  if (length(peps) > 1) {
    combos <- combn(peps, 2)
    peps <- unlist(lapply(seq_len(dim(combos)[2]), function(x) {
      return(remove_substr(combos[, x]))
    })) %>% unique()
  }
  align <- Biostrings::pairwiseAlignment(
    pattern = peps, subject = protein,
    type = "global-local"
  )
  vis <- sapply(aligned(align), as.character, USE.NAMES = FALSE)
  if (length(vis) > 1) {
    vis <- purrr::reduce(vis, resolveAlignment)
  }
  cov <- Biostrings::coverage(align)
  cov <- tibble(values = cov@values, lengths = cov@lengths)
  sum_cov <- cov %>%
    dplyr::filter(values > 0) %>%
    select(lengths) %>%
    sum()
  chars <- nchar(protein)
  nmatch <- sum(nmatch(align)) / chars
  if (nmatch >= 1) {
    # Redundant peptides (that differ by at least one residue) can cause nmatch to
    # exceed 1
    # We figure out the actual number of matches by taking the ratio of indices
    # that are not "-" in the alignment string
    total <- nchar(vis)
    nmatch <- (total - str_count(vis, "-")) / total
  }
  result <- list(
    tb = tibble(
      pcoverage_nmatch = nmatch,
      pcoverage_align = sum_cov / chars,
      alignment = vis
    ),
    seqs = aligned(align)
  )
  return(result)
}


ALIGNMENT_TB <- tibble()
appliedCoverage <- function(row) {
  peps <- str_split_1(row[["unique_peptides"]], ";")
  prot <- row[["seq"]]
  cov <- coverage(prot, peps)
  alignment_row <- tibble(
    alignment = purrr::map_chr(cov$seqs, as.character),
    ProteinId = row[["ProteinId"]],
    UniProtKB_ID = row[["UniProtKB_ID"]]
  )
  ALIGNMENT_TB <<- dplyr::bind_rows(ALIGNMENT_TB, alignment_row)
  return(cov$tb)
}

coverageCalc <- function(prot_df) {
  # Coverage calculation will only be performed on full-length proteins
  considered <- dplyr::filter(prot_df, !grepl("U|D", ProteinId))
  cov_info <- apply(considered, 1, appliedCoverage) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(ProteinId = considered$ProteinId)
  return(cov_info)
}

#' Split the dataframe to take advantage of parallelism
#'
splitForCoverage <- function(tb, path) {
  tb <- dplyr::filter(tb, !is.na(seq)) %>% dplyr::select(c(ProteinId, unique_peptides, seq, UniProtKB_ID))
  size <- 500
  num_rows <- nrow(tb)
  window <- c(0, size)
  slices <- list()
  while (num_rows >= size) {
    num_rows <- num_rows - size
    slices <- c(slices, list(dplyr::slice(tb, window[1]:window[2])))
    window <- window + size
  }
  slices <- c(slices, list(dplyr::slice(tb, window[1]:window[2])))
  for (n in seq_along(slices)) {
    write_tsv(slices[[n]], glue("{path}/combined_slice{n}.tsv"))
  }
}

mergeWithPattern <- function(path, pattern) {
  list.files(path, full.names = TRUE, pattern = pattern) %>%
    map(., read_tsv) %>%
    dplyr::bind_rows()
}

mergeCoverage <- function(tb_path, merge_into) {
  coverage <- mergeWithPattern(tb_path, ".*_calculated.tsv") %>%
    distinct()
  alignments <- mergeWithPattern(tb_path, ".*_alignments.tsv") %>%
    distinct()
  tb <- inner_join(merge_into, coverage, by = join_by(ProteinId))
  return(list(coverage = tb, alignment = alignments))
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"), type = "character")
  parser <- add_option(parser, c("-o", "--output_file"),
    type = "character",
    help = "Output file name"
  )
  parser <- add_option(parser, c("-p", "--output_path"),
    type = "character",
    help = "Output path"
  )
  parser <- add_option(parser, "--input_path")
  parser <- add_option(parser, c("-s", "--split"),
    action = "store_true", default = FALSE
  )
  parser <- add_option(parser, "--alignment_fasta",
    type = "character"
  )
  parser <- add_option(parser, "--alignment_tsv",
    type = "character"
  )
  parser <- add_option(parser, c("-c", "--calculate"),
    action = "store_true", default = FALSE
  )
  parser <- add_option(parser, c("-m", "--merge"),
    action = "store_true", default = FALSE
  )
  args <- parse_args(parser)
  input <- read_tsv(args$input)
  if (args$split) {
    splitForCoverage(input, args$output_path)
  } else if (args$calculate) {
    calculated <- coverageCalc(input)
    new_name <- gsub(
      "\\.tsv",
      "_calculated.tsv", args$input
    )
    alignment_file <- gsub(
      "\\.tsv",
      "_alignments.tsv", args$input
    )
    write_tsv(ALIGNMENT_TB, alignment_file)
    write_tsv(calculated, new_name)
  } else if (args$merge) {
    m <- mergeCoverage(args$input_path, input)
    apply(dplyr::filter(m$coverage, !is.na(seq)), 1, writeAlignments,
      file_name = args$alignment_fasta
    )
    merged <- dplyr::select(m$coverage, -alignment)
    write_tsv(merged, args$output_path)
    write_tsv(m$alignment, args$alignment_tsv)
  }
}
