library(tidyverse)
library(parallel)
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


resolveAlignment2 <- function(seq1, seq2) {
  if (length(seq1) > 1) {
    extra <- seq1[[2]]
  } else {
    extra <- NULL
  }
  if (class(seq1) == "list") {
    seq1 <- seq1[[1]]
  }
  algns <- list(c())
  split1 <- str_split_1(seq1, "")
  split2 <- str_split_1(seq2, "")
  sapply(seq_along(split1), \(x) {
    if (split1[x] == "-" && split2[x] == "-") {
      algns <<- appendAllList(algns, "-")
    } else if (split1[x] != "-" && split2[x] == "-") {
      algns <<- appendAllList(algns, split1[x])
    } else if (split1[x] == "-" && split2[x] != "-") {
      algns <<- appendAllList(algns, split2[x])
    } else if (split1[x] == split2[x]) {
      algns <<- appendAllList(algns, split2[x])
    } else {
      if (length(algns) == 1) {
        algns[[2]] <<- algns[[1]]
      }
      algns <<- appendAllList(algns, c(split1[x], split2[x]))
    }
  })
  algns <- c(algns, extra)
  algns <- lapply(algns, paste0, collapse = "")
  return(algns)
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
    filter(values > 0) %>%
    select(lengths) %>%
    sum()
  chars <- nchar(protein)
  nmatch <- sum(nmatch(align)) / chars
  result <- tibble(
    coverage_nmatch = ifelse(nmatch >= 1, 1.00, nmatch),
    coverage_alignlen = sum_cov / chars,
    alignment = vis
  )
  return(result)
}

appliedCoverage <- function(row) {
  peps <- str_split_1(row[["unique_peptides"]], ",")
  if (length(peps) > 1) {
    pairs <- combn(peps, 2)
    peps <- unlist(lapply(seq_len(dim(pairs)[2]), function(x) {
      remove_substr(pairs[, x])
    }), use.names = FALSE)
  }
  prot <- row[["seq"]]
  return(coverage(prot, peps))
}


coverageCalc <- function(prot_df, parallel) {
  # Coverage calculation will only be performed on full-length proteins
  not_considered <- filter(prot_df, is.na(seq) | grepl("U|D", ProteinId))
  considered <- filter(prot_df, !is.na(seq) & !grepl("U|D", ProteinId))
  cov_info <- apply(considered, 1, appliedCoverage) %>% bind_rows()
  considered <- bind_cols(considered, cov_info)
  final <- bind_rows(considered, not_considered)
  return(final)
}
