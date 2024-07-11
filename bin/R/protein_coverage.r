library(tidyverse)
library(seqinr)
library(glue)
library(Biostrings)


write_alignments <- function(row, file_name) {
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

split_with_groups <- function(str) {
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
resolve_alignment <- function(seq1, seq2) {
  if (grepl("\\[", seq1)) {
    split1 <- split_with_groups(seq1)
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

append_all_list <- function(lst, to_append) {
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

get_coverage <- function(protein, peps) {
  # Using alignment percent identity for the coverage
  # calculation takes into account gaps and mismatches
  ids <- names(peps)
  align <- Biostrings::pairwiseAlignment(
    pattern = peps, subject = protein,
    type = "global-local"
  )
  vis <- sapply(aligned(align), as.character, USE.NAMES = FALSE)
  if (length(vis) > 1) {
    vis <- purrr::reduce(vis, resolve_alignment)
  }
  cov <- Biostrings::coverage(align)
  # "Lengths" in the Rle object are stretches of peptide residues
  # "Value" indicates whether or not it was covered
  cov <- tibble(values = cov@values, lengths = cov@lengths)
  sum_cov <- cov %>%
    dplyr::filter(values > 0) %>%
    select(lengths) %>%
    sum()
  chars <- nchar(protein)
  result <- list(
    tb = tibble(
      pcoverage_align = sum_cov / chars,
      alignment = vis,
    ),
    seqs = tibble(
      alignment = purrr::map_chr(aligned(align), as.character),
      original = peps,
      id = ids,
      start = start(subject(align)),
      end = end(subject(align)),
      score = score(align),
      n_match = nmatch(align),
      n_mismatch = nmismatch(align),
      seq_length = chars
    )
  )
  return(result)
}


get_id_from_rows <- function(row, id_col) {
  if (nrow(row) > 1) {
    unmatched_removed <- row |> filter(!grepl("D", !!as.symbol(id_col)))
    if (nrow(unmatched_removed) > 0) {
      row <- slice_sample(unmatched_removed, n = 1)
    } else {
      row <- slice_sample(row, n = 1)
    }
  }
  return(row[[id_col]])
}



# `seq_map` is `seq-header_map_found.tsv`
# `peptide_map` is `percolator_peptide_map.tsv`
# both from combine_percolator
identify_peptides <- function(protein_id, matched_peptide_ids, unique_peptides) {
  if (is.na(matched_peptide_ids)) {
    matched <- NA
  } else {
    matched <- str_split_1(matched_peptide_ids, ";")
  }
  peps <- unique_peptides |> str_split_1(";")
  cur_seq_map <- SEQ_MAP |> filter(id %in% matched)
  cur_peptide_map <- PEPTIDE_MAP |> filter(ProteinId == protein_id | ProteinId %in% matched)
  pep_metadata <- peps |>
    map_chr(
      \(x) {
        # Want to search for the id in a narrowed window so as to reduce instances
        # of spurious matches
        id_row <- filter(cur_seq_map, seq == x)
        if (nrow(id_row) > 0) {
          return(get_id_from_rows(id_row, "id"))
        }
        id_row <- filter(cur_peptide_map, peptideIds == x)
        if (nrow(id_row) > 0) {
          return(get_id_from_rows(id_row, "ProteinId"))
        }
        id_row <- filter(PEPTIDE_MAP, peptideIds == x)
        if (nrow(id_row) > 0) {
          return(get_id_from_rows(id_row, "ProteinId"))
        }
        id_row <- filter(SEQ_MAP, seq == x)
        if (nrow(id_row) > 0) {
          return(get_id_from_rows(id_row, "id"))
        }
        "unknown"
      }
    )
  names(peps) <- pep_metadata
  peps
}


ALIGNMENT_TB <- tibble()
applied_coverage <- function(
    protein_id, matched_peptide_ids, unique_peptides, seq,
    uniprot_id) {
  peps <- str_split_1(unique_peptides, ";")
  peps <- identify_peptides(protein_id, matched_peptide_ids, unique_peptides)
  cov <- get_coverage(seq, peps)
  alignment_row <- cov$seqs |> mutate(
    ProteinId = protein_id,
    UniProtKB_ID = uniprot_id
  )
  ALIGNMENT_TB <<- dplyr::bind_rows(ALIGNMENT_TB, alignment_row)
  return(cov$tb)
}

calculate_coverage <- function(prot_df) {
  # Coverage calculation will only be performed on full-length proteins
  considered <- dplyr::filter(prot_df, !grepl("U|D", ProteinId))
  cov_info <- pmap(
    list(
      considered$ProteinId,
      considered$MatchedPeptideIds,
      considered$unique_peptides,
      considered$seq,
      considered$UniProtKB_ID
    ),
    applied_coverage
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(ProteinId = considered$ProteinId)
  return(cov_info)
}

#' Split the dataframe to take advantage of parallelism
#'
split_for_coverage <- function(tb, path) {
  tb <- dplyr::filter(tb, !is.na(seq)) %>% dplyr::select(c(
    ProteinId,
    unique_peptides, seq, UniProtKB_ID, MatchedPeptideIds
  ))
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

merge_with_pattern <- function(path, pattern) {
  list.files(path, full.names = TRUE, pattern = pattern) %>%
    map(., read_tsv) %>%
    dplyr::bind_rows()
}

#' Merge coverage results back into overall results
#'
#' @param merge_into tb containing overall results
#' @param tb_path path to coverage files
merge_coverage <- function(tb_path, merge_into) {
  coverage <- merge_with_pattern(tb_path, ".*_calculated.tsv") %>%
    distinct()
  alignments <- merge_with_pattern(tb_path, ".*_alignments.tsv") %>%
    distinct()
  tb <- left_join(merge_into, coverage, by = join_by(ProteinId))
  # Left join so as not to lose identified de novo peptides (which don't have coverage)
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
  parser <- add_option(parser, "--seq_header_map")
  parser <- add_option(parser, "--peptide_map")
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
    split_for_coverage(input, args$output_path)
  } else if (args$calculate) {
    SEQ_MAP <- read_tsv(args$seq_header_map)
    PEPTIDE_MAP <- read_tsv(args$peptide_map)
    calculated <- calculate_coverage(input)
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
    m <- merge_coverage(args$input_path, input)
    apply(dplyr::filter(m$coverage, !is.na(seq)), 1, write_alignments,
      file_name = args$alignment_fasta
    )
    merged <- dplyr::select(m$coverage, -alignment)
    write_tsv(merged, args$output_path)
    write_tsv(m$alignment, args$alignment_tsv)
  }
}
