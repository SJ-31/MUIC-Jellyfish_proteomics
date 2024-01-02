library(tidyverse)
library(Peptides)
library(optparse)
library(glue)

read_percolator <- function(file) {
  t <- read_tsv(file, col_names = FALSE) %>%
    mutate(X6 = unlist(lapply(X6, gsub, pattern = "\t", replacement = ";"))) %>%
    `colnames<-`(.[1, ]) %>%
    slice(-1) %>%
    select(c(PSMId, peptide, proteinIds)) %>%
    rename(protein = proteinIds)
  return(t)
}

split_ambiguous <- function(table, index) {
  base <- strsplit(table[index, ]$`Base.Sequence`, "\\|")[[1]]
  row <- lapply(seq_along(base), function(x) {
    return(data.frame(
      `File.Name` = table[index, ]$`File.Name`,
      `Precursor.Charge` = table[index, ]$`Precursor.Charge`,
      `Base.Sequence` = base[x],
      `Protein.Accession` = table$`Protein.Accession`,
      `Scan.Number` = table$`Scan.Number`
    ))
  })[[1]]
  return(row)
}

sort_ambiguous <- function(mm) {
  #' Split a metamorpheus psm matching to multiple proteins into multiple rows
  ambiguous <- mm %>% filter(grepl("\\|", `Base.Sequence`))
  if (dim(ambiguous)[1] == 0) {
    return(mm)
  }
  single <- mm %>% filter(!grepl("\\|", `Base.Sequence`))
  changed <- lapply(1:dim(ambiguous)[1], split_ambiguous,
    table = ambiguous
  ) %>%
    bind_rows()
  bound <- bind_rows(list(single, changed)) %>%
    as_tibble()
  return(bound)
}

read_metamorpheus <- function(metamorpheus_file) {
  # Needs metamorpheus AllPSMS.psmtsv
  old_names <- c(
    "File.Name", "Precursor.Charge",
    "Base.Sequence",
    "Protein.Accession",
    "Scan.Number"
  )
  new_names <- c("file", "precursorCharge", "peptide", "protein", "scan")
  mm <- read.delim(metamorpheus_file, sep = "\t") %>%
    as_tibble() %>%
    filter(Decoy == "N") %>%
    select(all_of(old_names))
  mm <- sort_ambiguous(mm) %>%
    distinct(`Base.Sequence`, .keep_all = TRUE) %>%
    mutate(`Scan.Number` = paste0(
      `File.Name`, ".",
      `Scan.Number`
    )) %>%
    rename_with(~new_names, all_of(old_names)) %>%
    mutate(protein = unlist(lapply(protein, gsub,
      pattern = "\\|",
      replacement = ";"
    )))
  return(mm)
}


read_tide <- function(tide_file, mapping) {
  # Needs tide-search.target.txt
  selection <- c("file", "scan", "charge", "sequence", "protein.id")
  tide <- read.delim(tide_file, sep = "\t")
  t <- tide %>%
    select(all_of(selection)) %>%
    mutate(file = unlist(lapply(file,
      gsub,
      pattern = "\\..*",
      replacement = ""
    ))) %>%
    mutate(scan = unlist(lapply(seq_along(scan), function(x) {
      return(paste0(.[x, ]$file, ".", .[x, ]$scan))
    }))) %>%
    mutate(protein = unlist(lapply(
      seq_along(scan),
      function(x) {
        grouped <- .[x, ]$`protein.id` %>%
          str_split(",", simplify = TRUE) %>%
          lapply(., gsub, pattern = "\\(.*\\)", replacement = "") %>%
          paste0(collapse = ";")
        return(grouped)
      }
    ))) %>%
    rename(peptide = sequence) %>%
    inner_join(., mapping, by = join_by(scan == scanNum)) %>%
    as_tibble()
  return(t)
}

##  Functions for formatting scan number
##    The scan number will let you match retention time
comet_scans <- function(comet_id_str) {
  path_removed <- gsub(".*/", "", comet_id_str)
  name_cleared <- strsplit(path_removed, "_") %>% unlist(use.names = FALSE)
  return(glue("{name_cleared[1]}.{name_cleared[2]}"))
}

msfragger_scans <- function(msfragger_id_str) {
  match_result <- regexec("(.*\\..*)\\..*\\..*", msfragger_id_str)
  groups <- regmatches(msfragger_id_str, match_result)[[1]]
  return(groups[2])
}

msgf_scans <- function(msgf_id_str) {
  return(gsub("_SII", ".", msgf_id_str))
}

get_file_name <- function(scan) {
  return(gsub("\\..*", "", scan))
}

clean_peptide <- function(modified_pep) {
  clean_pep <- str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
    paste0(collapse = "")
  return(clean_pep)
}

read_engine_psms <- function(args) {
  percolator_input <- args$input
  engine <- args$engine
  mapping <- read.delim(args$msms_mapping, sep = "\t")
  if (grepl("metamorpheus", engine)) {
    psms <- read_metamorpheus(percolator_input)
  } else if (engine == "tide") {
    psms <- read_tide(percolator_input, mapping)
  } else {
    psms <- read_percolator(percolator_input)
    if (engine == "comet") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, comet_scans)))
    } else if (engine == "msgf") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, msgf_scans)))
    } else if (grepl("msfragger", engine)) {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, msfragger_scans)))
    } else if (engine == "identipy") {
      psms <- psms %>% rename(scan = PSMId)
    }
  }
  psms <- psms %>%
    mutate(file = unlist(lapply(scan, get_file_name))) %>%
    mutate(base_peptide = unlist(lapply(peptide, clean_peptide))) %>%
    select(c(
      "file", "scan", "base_peptide", "peptide",
      "protein"
    )) %>%
    mutate(mw = unlist(lapply(base_peptide, mw, monoisotopic = TRUE)))
  psms <- inner_join(psms, mapping, by = join_by(scan == scanNum)) %>%
    mutate(engine = engine) %>%
    mutate(protein = unlist(lapply(protein, gsub, pattern = '"', replacement = "")))
  return(psms)
}

expandProteinRows <- function(row) {
  ProteinId <- row[["ProteinId"]] %>% str_split_1(",")
  peptide <- row[["peptideIds"]] %>% str_split_1(" ")
  combos <- crossing(ProteinId, peptide)
  return(combos)
}

noTermini <- function(peptide_col) {
  unlist(lapply(peptide_col, \(x) {
    return(substr(x, start = 3, stop = nchar(x) - 2))
  }))
}

merge_unmatched <- function(final_df, unmatched_peptides, proteins) {
  #' Merge the final scan numbers with the ids of the unmatched peptides
  #' so that they may be quantified
  prot_df <- read_tsv(proteins) %>%
    select(ProteinId, peptideIds) %>%
    apply(., 1, expandProteinRows) %>%
    bind_rows()
  u_df <- read_tsv(unmatched_peptides)
  unwanted_cols <- colnames(u_df) %>% discard(., \(x) x %in% "peptideIds")
  # Join with unmatched peptides
  no_prot1 <- final_df %>%
    filter(is.na(protein)) %>%
    mutate(no_termini = noTermini(peptide))
  has_prot1 <- final_df %>% filter(!is.na(protein))
  joined1 <- left_join(no_prot1, u_df,
    by =
      join_by(x$no_termini == y$peptideIds)
  ) %>%
    mutate(protein = ProteinId) %>%
    select(-all_of(unwanted_cols))
  # Join with other proteins matched by percolator
  no_prot2 <- joined1 %>%
    filter(is.na(protein))
  has_prot2 <- joined1 %>% filter(!is.na(protein))
  joined2 <- left_join(no_prot2, prot_df,
    by =
      join_by(x$no_termini == y$peptide)
  ) %>%
    mutate(protein = ProteinId) %>%
    select(-ProteinId)
  return(select(
    bind_rows(joined1, has_prot1, joined2, has_prot2),
    -no_termini
  ))
}


if (sys.nframe() == 0) {
  parser <- OptionParser()
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  parser <- add_option(parser, c("-m", "--msms_mapping"),
    type = "character",
    help = "MsMs mapping file"
  )
  parser <- add_option(parser, c("-i", "--input"),
    type = "character",
    help = "Percolator psm file"
  )
  parser <- add_option(parser, c("-p", "--protein"),
    type = "character",
    help = "Percolator protein file"
  )
  parser <- add_option(parser, c("-u", "--unmatched_peptides"),
    type = "character",
    help = "TSV containing unmatched peptides"
  )
  parser <- add_option(parser, c("-e", "--engine"),
    type = "character",
    help = "Engine psm file was obtained from"
  )
  args <- parse_args(parser)
  f <- read_engine_psms(args)
  if (any(is.na(f$protein))) {
    final <- merge_unmatched(f, args$unmatched_peptides, args$protein)
    write_delim(final, args$output, delim = "\t")
  } else {
    write_delim(f, args$output, delim = "\t")
  }
}
