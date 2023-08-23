library(tidyverse)
library(Peptides)
library(optparse)
library(glue)
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)

flashlfq_header <- c(
  "File Name", "Scan Retention Time", "Precursor Charge",
  "Base Sequence", "Full Sequence", "Peptide Monoisotopic Mass",
  "Protein Accession"
)

clean_peptide <- function(modified_pep) {
  return(str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
    paste0(collapse = ""))
}

group_prot <- function(prot_vec) {
  grouped <- paste0(prot_vec, collapse = ";")
  split_groups <- grouped %>%
    strsplit(";") %>%
    unlist(use.names = FALSE) %>%
    unique()
  return(paste0(split_groups, collapse = "|"))
}

merge_with_msms <- function(table, mapping) {
  return(table %>%
         mutate(`Peptide Monoisotopic Mass` = unlist(
      lapply(`Base Sequence`, mw, monoisotopic = TRUE)
    )) %>%
    inner_join(., mapping, by = join_by(scan == scanNum)) %>%
    rename("Scan Retention Time" = "retensionTime") %>%
    rename("Precursor Charge" = "precursorCharge") %>%
    select(all_of(flashlfq_header)) %>%
    as_tibble())
}

get_percolator_row <- function(row_index, percolator_lines) {
  splits <- percolator_lines[row_index] %>%
    str_split("\t") %>%
    unlist()
  return(tibble(
    PSMId = splits[1],
    peptide = splits[5],
    `Protein Accession` = group_prot(splits[6:length(splits)])
  ))
}

read_percolator <- function(percolator_file) {
  lines <- read_lines(percolator_file)
  p_tibble <- lapply(seq_along(lines)[-1], get_percolator_row,
    percolator_lines = lines
  ) %>% bind_rows()
  return(p_tibble)
}

get_maxquant_row <- function(row_index, maxquant_lines) {
  splits <- maxquant_lines[row_index] %>%
    str_split("\t") %>%
    unlist()
  return(tibble(
    `File Name` = str_match(splits[1], "(.*)\\..*\\..*")[, 2],
    label = splits[2],
    scan = splits[3],
    `Full Sequence` = splits[104],
    `Base Sequence` = clean_peptide(splits[104]),
    `Protein Accession` = group_prot(splits[105:length(splits)])
  ))
}

read_maxquant <- function(maxquant_file, mapping) {
  # Needs maxquant_all_pins
  lines <- read_lines(maxquant_file)
  q_tibble <- lapply(seq_along(lines)[-1], get_maxquant_row,
    maxquant_lines = lines
  ) %>%
    bind_rows() %>%
    filter(label == 1) %>%
    filter(!(grepl("REV__", `Protein Accession`, ))) %>%
    select(-label) %>%
    mutate(scan = paste0(`File Name`, ".", scan)) %>%
    merge_with_msms(mapping)
  return(q_tibble)
}

split_ambiguous <- function(table, index) {
  base <- strsplit(table[index, ]$`Base.Sequence`, "\\|")[[1]]
  mods <- strsplit(table[index, ]$`Full.Sequence`, "\\|")[[1]]
  mass <- table[index, ]$`Peptide.Monoisotopic.Mass` %>% gsub("\\|.*", "", .)
  row <- lapply(seq_along(base), function(x) {
    return(data.frame(
      `File.Name` = table[index, ]$`File.Name`,
      `Scan.Retention Time` = table[index, ]$`Scan.Retention.Time`,
      `Precursor.Charge` = table[index, ]$`Precursor.Charge`,
      `Base.Sequence` = base[x],
      `Full.Sequence` = mods[x],
      `Peptide.Monoisotopic.Mass` = as.double(mass),
      `Protein.Accession` = table$`Protein.Accession`
    ))
  })[[1]]
  return(row)
}

sort_ambiguous <- function(mm) {
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
  # Needs metamorpheus AllPSMs.psmtsv
  old_names <- c(
    "File.Name", "Scan.Retention.Time", "Precursor.Charge",
    "Base.Sequence", "Full.Sequence", "Peptide.Monoisotopic.Mass",
    "Protein.Accession"
  )
  new_names <- lapply(old_names, gsub,
    pattern = "\\.",
    replacement = " "
  ) %>% unlist(use.names = FALSE)
  mm <- read.delim(metamorpheus_file, sep = "\t") %>%
    as_tibble() %>%
    filter(Decoy == "N") %>%
    select(all_of(old_names)) %>%
    mutate(
      `Peptide.Monoisotopic.Mass` =
        as.double(`Peptide.Monoisotopic.Mass`)
    )
  mm <- sort_ambiguous(mm) %>%
    rename_with(~new_names, all_of(old_names)) %>%
    distinct(`Base Sequence`, .keep_all = TRUE)
  return(mm)
}


read_tide <- function(tide_file, mapping) {
  # Needs tide-search.target.txt
  old_cols <- c("file", "scan", "sequence", "protein.id")
  new_cols <- c(
    "File Name", "scan",
    "Full Sequence", "Protein Accession"
  )
  tide <- read.delim(tide_file, sep = "\t")
  t <- tide %>%
    select(all_of(old_cols)) %>%
    rename_with(~new_cols, all_of(old_cols)) %>%
    mutate(`File Name` = unlist(lapply(`File Name`,
      gsub,
      pattern = "\\..*",
      replacement = ""
    ))) %>%
    mutate(scan = unlist(lapply(seq_along(scan), function(x) {
      return(paste0(.[x, ]$`File Name`, ".", .[x, ]$scan))
    }))) %>%
    mutate(`Protein Accession` = unlist(lapply(
      seq_along(scan),
      function(x) {
        grouped <- .[x, ]$`Protein Accession` %>%
          str_split(",", simplify = TRUE) %>%
          lapply(., gsub, pattern = "\\(.*\\)", replacement = "") %>%
          paste0(collapse = "|")
        return(grouped)
      }
    ))) %>%
    mutate(`Base Sequence` = unlist(lapply(
      `Full Sequence`,
      clean_peptide
    ))) %>%
  merge_with_msms(mapping = mapping)
  return(t)
}


##  Functions for formatting scan number
##    The scan number will let you match retension time
comet_scans <- function(comet_id_str) {
  path_removed <- gsub(".*/", "", comet_id_str)
  name_cleared <- strsplit(path_removed, "_") %>% unlist(use.names = FALSE)
  return(glue("{name_cleared[1]}.{name_cleared[2]}"))
}

maxquant_scans <- function(mq_id_str) {
  match_result <- regexec("(.*\\..*)\\..*", mq_id_str)
  groups <- regmatches(mq_id_str, match_result)[[1]]
  return(groups[2])
}

msfragger_scans <- function(msfragger_id_str) {
  match_result <- regexec("(.*\\..*)\\..*\\..*", msfragger_id_str)
  groups <- regmatches(msfragger_id_str, match_result)[[1]]
  return(groups[2])
}

get_file_name <- function(scan) {
  return(gsub("\\..*", "", scan))
}


read_engine_psms <- function(percolator_input, engine, mapping) {
  if (engine == "metamorpheus") {
    return(read_metamorpheus(percolator_input))
  } else if (engine == "tide") {
    return(read_tide(percolator_input, mapping))
  } else if (engine == "maxquant") {
    return(read_maxquant(percolator_input, mapping))
  }
  # Read the percolator psm file from <engine> and format the
  #     the output according to flashlfq
  #     This relies on all functions defined above
  psms <- read_percolator(percolator_input)
  if (engine == "comet") {
    psms <- psms %>% mutate(scan = unlist(lapply(PSMId, comet_scans)))
  } else if (engine == "maxquant") {
    psms <- psms %>% mutate(scan = unlist(lapply(PSMId, maxquant_scans)))
  } else if (engine == "msfragger") {
    psms <- psms %>% mutate(scan = unlist(lapply(PSMId, msfragger_scans)))
  } else if (engine == "identipy") {
    psms <- psms %>% mutate(scan = PSMId)
  }
  psms <- psms %>%
    mutate(`File Name` = unlist(lapply(scan, get_file_name))) %>%
    mutate(`Full Sequence` = peptide) %>%
    mutate(`Base Sequence` = unlist(lapply(peptide, clean_peptide))) %>%
    mutate(`Peptide Monoisotopic Mass` = unlist(
      lapply(`Base Sequence`, mw, monoisotopic = TRUE)
    )) %>%
    select(c(
      "File Name", "scan", "Base Sequence", "Full Sequence",
      "Peptide Monoisotopic Mass",
      "Protein Accession"
    ))
  new_cols <- c("Scan Retention Time", "Precursor Charge")
  old_cols <- c("retensionTime", "precursorCharge")
  psms <- inner_join(psms, mapping, by = join_by(scan == scanNum)) %>%
    rename_with(~new_cols, all_of(old_cols)) %>%
    select(all_of(flashlfq_header))
  return(psms)
}

parser <- OptionParser()
parser <- add_option(parser, c("--metamorpheus"),
  type = "character",
  help = "metamorpheus psms"
)
parser <- add_option(parser, c("--identipy"),
  type = "character",
  help = "identipy psms"
)
parser <- add_option(parser, c("--maxquant"),
  type = "character",
  help = "maxquant psms"
)
parser <- add_option(parser, c("--msfragger"),
  type = "character",
  help = "msfragger psms"
)
parser <- add_option(parser, c("--comet"),
  type = "character",
  help = "comet psms"
)
parser <- add_option(parser, c("--tide"),
  type = "character",
  help = "tide psms"
)
parser <- add_option(parser, c("-o", "--output"),
  type = "character",
  help = "Output file name"
)
parser <- add_option(parser, c("-m", "--msms_mapping"),
  type = "character",
  help = "MsMs mapping file"
)
args <- parse_args(parser)

mapping <- read.delim(args$msms_mapping, sep = "\t")


file_list <- list(
  comet = args$comet, identipy = args$identipy,
  tide = args$tide, metamorpheus = args$metamorpheus,
  msfragger = args$msfragger, maxquant = args$maxquant
)

# for testing
## file_list <- list(
##   comet = "../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv",
##   identipy = "../results/test_manifest/1-First_pass/Percolator/identipy_percolator_psms.tsv", maxquant = "../results/test_manifest/1-First_pass/Percolator/maxquant_percolator_psms.tsv", msfragger = "../results/test_manifest/1-First_pass/Percolator/msfragger_percolator_psms.tsv",
##   metamorpheus = "../results/test_manifest/1-First_pass/Metamorpheus/metamorpheus_AllPSMs.psmtsv", tide = "../results/test_manifest/1-First_pass/Tide/tide-search.target.txt"
## )
## mapping <- read.delim("../results/CiCs_metrics.tsv", sep = "\t")


flashlfq <- read_engine_psms(file_list[[1]], names(file_list)[1], mapping)
# Set up the first set of psms
unique_peps <- flashlfq$`Base Sequence`
#
# Add only new psms to the first set
for (n in seq_along(file_list)[-1]) {
  current <- read_engine_psms(file_list[[n]], names(file_list)[n], mapping)
  current_peps <- current$`Base Sequence`
  new_rows <- current %>% filter(!(`Base Sequence` %in% unique_peps))
  flashlfq <- bind_rows(flashlfq, new_rows)
  unique_peps <- c(
    unique_peps,
    current_peps[!(current_peps %in% unique_peps)]
  )
  rm(current, current_peps, new_rows)
}

flashlfq <- flashlfq %>%
  filter(!(is.na(`Peptide Monoisotopic Mass`))) %>%
  mutate(`Base Sequence` = unlist(lapply(`Base Sequence`, gsub, pattern = "X", replacement = "")))

write_delim(flashlfq, args$output, delim = "\t")
