library(MSnbase)
library(tidyverse)
library(optparse)
library(glue)

group_prot <- function(prot_vec) {
  grouped <- paste0(prot_vec, collapse = ";")
  split_groups <- grouped %>%
    strsplit(";") %>%
    unlist(use.names = FALSE) %>%
    unique()
  return(paste0(split_groups, collapse = ";"))
}

get_percolator_row <- function(row_index, percolator_lines) {
  splits <- percolator_lines[row_index] %>%
    str_split("\t") %>%
    unlist()
  return(tibble(
    PSMId = splits[1],
    peptide = splits[5],
    protein = group_prot(splits[6:length(splits)])
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
    file = str_match(splits[1], "(.*)\\..*\\..*")[, 2],
    label = splits[2],
    scan = splits[3],
    peptide = splits[104],
    protein = group_prot(splits[105:length(splits)])
  ))
}

read_maxquant <- function(maxquant_file, mapping) {
  lines <- read_lines(maxquant_file)
  q_tibble <- lapply(seq_along(lines)[-1], get_maxquant_row,
    maxquant_lines = lines
  ) %>%
    bind_rows() %>%
    filter(label == 1) %>%
    filter(!(grepl("REV__", protein, ))) %>%
    select(-label) %>%
    mutate(scan = paste0(file, ".", scan))
  return(q_tibble)
}

read_metamorpheus <- function(metamorpheus_file, mapping) {
  old_names <- c(
    "File.Name",
    "Precursor.Scan.Number"
  )
  new_names <- c("file", "scan")
  mm <- read.delim(metamorpheus_file, sep = "\t") %>%
    as_tibble() %>%
    filter(Decoy == "N") %>%
    select(all_of(old_names)) %>%
    mutate(`Precursor.Scan.Number` = paste0(
      `File.Name`, ".",
      `Precursor.Scan.Number`
    )) %>%
    rename_with(~new_names, all_of(old_names))
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
    dplyr::rename(peptide = sequence)
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


read_engine_psms <- function(percolator_input, engine, mapping) {
  if (engine == "metamorpheus") {
    psms <- read_metamorpheus(percolator_input, mapping)
  } else if (engine == "tide") {
    psms <- read_tide(percolator_input, mapping)
  } else if (engine == "maxquant") {
    psms <- read_maxquant(percolator_input, mapping)
  } else {
    psms <- read_percolator(percolator_input)
    if (engine == "comet") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, comet_scans)))
    } else if (engine == "maxquant") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, maxquant_scans)))
    } else if (engine == "msfragger") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, msfragger_scans)))
    } else if (engine == "identipy") {
      psms <- psms %>% dplyr::rename(scan = PSMId)
    }
  }
  psms <- psms %>%
    inner_join(., mapping, by = join_by(scan == scanNum)) %>%
    select(c("scan")) %>%
    unique()
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
parser <- add_option(parser, c("-m", "--msms_mapping"),
  type = "character",
  help = "MsMs mapping file"
)
parser <- add_option(parser, c("--mzML_path"),
  type = "character",
  help = "Path to searched mzML files"
)
args <- parse_args(parser)
mapping <- read.delim(args$msms_mapping, sep = "\t")

file_list <- list(
  comet = args$comet, identipy = args$identipy,
  tide = args$tide, metamorpheus = args$metamorpheus,
  msfragger = args$msfragger, maxquant = args$maxquant
)

## file_list <- list( # Testing
##   comet = "../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv",
##   identipy = "../results/test_manifest/1-First_pass/Percolator/identipy_percolator_psms.tsv", maxquant = "../results/test_manifest/1-First_pass/MaxQuant/maxquant_all_pins.temp", msfragger = "../results/test_manifest/1-First_pass/Percolator/msfragger_percolator_psms.tsv",
##   metamorpheus = "../results/test_manifest/1-First_pass/Metamorpheus/metamorpheus_AllPSMs.psmtsv", tide = "../results/test_manifest/1-First_pass/Tide/tide-search.target.txt"
## )


mapping <- read.delim("../results/CiCs_metrics.tsv", sep = "\t")


all_engines <- lapply(seq_along(file_list), function(x) {
  return(read_engine_psms(file_list[[x]], names(file_list)[x], mapping))
})
matched_psms <- all_engines %>% unlist() %>% unique()
rm(all_engines)

unmatched_msms <- mapping %>% filter(!(scanNum %in% matched_psms)) %>%
  select(scanNum)

unmatched_df <- seq_along(unmatched_msms[[1]]) %>% lapply(., function(x) {
    split <- strsplit(unmatched_msms[x, ], "\\.")[[1]]
    print(split)
    return(data.frame(file = split[1], scan = split[2]))
  }) %>%
  bind_rows() %>%
  mutate(scan = lapply(scan, function(x) {return(glue("F1.S{x}"))}))

path <- args$mzML_path
## path <- "../results" # For tests

msms_files <- as.list(list.files(path,
                                 pattern = ".*\\.mzML" )) %>%
  `names<-`(gsub("\\..*", "", unlist(.)))

filter_unmatched_msms <- function(file_name, msms_path, df) {
  df <- filter(df, file == file_name)
  mzml <- readMSData(msms)
  sp <- spectra(mzml)
  unmatched_mzml <- mzml[names(sp) %in% file_scans$scan]
  rm(mzml)
  return(unmatched_mzml)
}

unmatched <- lapply(names(msms_files), function(x) {
  filter_unmatched_msms(x, msms_files[[x]], unmatched_df)
}) %>% `names<-`(names(msms_files))

lapply(names(unmatched), function (x) {
  writeMSData(unmatched[[x]], file = glue("unmatched_{x}.mzML"))
})
