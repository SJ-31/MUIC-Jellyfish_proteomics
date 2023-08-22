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

read_maxquant <- function(maxquant_file) {
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

split_ambiguous <- function(table, index) {
  base <- strsplit(table[index, ]$`Base.Sequence`, "\\|")[[1]]
  row <- lapply(seq_along(base), function(x) {
    return(data.frame(
      `File.Name` = table[index, ]$`File.Name`,
      `Precursor.Charge` = table[index, ]$`Precursor.Charge`,
      `Base.Sequence` = base[x],
      `Protein.Accession` = table$`Protein.Accession`,
      `Precursor.Scan.Number` = table$`Precursor.Scan.Number`
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
  old_names <- c(
    "File.Name", "Precursor.Charge",
    "Base.Sequence",
    "Protein.Accession",
    "Precursor.Scan.Number"
  )
  new_names <- c("file", "precursorCharge", "peptide", "protein", "scan")
  mm <- read.delim(metamorpheus_file, sep = "\t") %>%
    as_tibble() %>%
    filter(Decoy == "N") %>%
    select(all_of(old_names))
  mm <- sort_ambiguous(mm) %>%
    distinct(`Base.Sequence`, .keep_all = TRUE) %>%
    mutate(`Precursor.Scan.Number` = paste0(
      `File.Name`, ".",
      `Precursor.Scan.Number`
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

clean_peptide <- function(modified_pep) {
  clean_pep <- str_extract_all(modified_pep, "[A-Z]+")[[1]] %>%
    paste0(collapse = "")
  return(glue("_{clean_pep}_"))
}

file_pivot <- function(psm_df) {
  return(psm_df %>%
    select(c("file", "protein", "ion", "precursorIntensity")) %>%
    distinct(ion, .keep_all = TRUE) %>%
    pivot_wider(names_from = file, values_from = precursorIntensity))
}

read_engine_psms <- function(percolator_input, engine, mapping) {
  if (engine == "metamorpheus") {
    psms <- read_metamorpheus(percolator_input)
  } else if (engine == "tide") {
    psms <- read_tide(percolator_input, mapping)
  } else if (engine == "maxquant") {
    psms <- read_maxquant(percolator_input)
  } else {
    psms <- read_percolator(percolator_input)
    if (engine == "comet") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, comet_scans)))
    } else if (engine == "maxquant") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, maxquant_scans)))
      # id != scan number fix it
    } else if (engine == "msfragger") {
      psms <- psms %>% mutate(scan = unlist(lapply(PSMId, msfragger_scans)))
    } else if (engine == "identipy") {
      psms <- psms %>% rename(scan = PSMId)
    }
  }
  psms <- psms %>%
    mutate(file = unlist(lapply(scan, get_file_name))) %>%
    mutate(ion = unlist(lapply(peptide, clean_peptide))) %>%
    select(c(
      "file", "scan", "ion",
      "protein"
    ))
  psms <- inner_join(psms, mapping, by = join_by(scan == scanNum)) %>%
    mutate(ion = paste0(ion, precursorCharge)) %>%
    file_pivot()
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

## file_list <- list(
##   comet = "../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv",
##   identipy = "../results/test_manifest/1-First_pass/Percolator/identipy_percolator_psms.tsv", maxquant = "../results/test_manifest/1-First_pass/MaxQuant/maxquant_all_pins.temp", msfragger = "../results/test_manifest/1-First_pass/Percolator/msfragger_percolator_psms.tsv",
##   metamorpheus = "../results/test_manifest/1-First_pass/Metamorpheus/metamorpheus_AllPSMs.psmtsv", tide = "../results/test_manifest/1-First_pass/Tide/tide-search.target.txt"
## )
## mapping <- read.delim("../results/CiCs_metrics.tsv", sep = "\t")

file_names <- mapping$scanNum %>%
  lapply(., gsub, pattern = "\\..*", replacement = "") %>%
  unlist(use.names = FALSE) %>%
  unique()

all_engines <- lapply(seq_along(file_list), function(x) {
  return(read_engine_psms(file_list[[x]], names(file_list)[x], mapping))
})
rm(all_engines)

# Merge matches between engines
joined <- Reduce(function(x, y) { full_join(x, y, by = "ion") }, all_engines)

# Group up proteins that share an ion
prot <- joined %>% select(grep("protein", colnames(joined)))
prot_col <- lapply(seq_along(1:dim(prot)[1]), function(x) {
  filtered <- prot[x, ][grepl("[a-zA-Z1-9]+", prot[1, ])]
  return(paste0(filtered, collapse = ";"))
}) %>%
  unlist()
joined <- joined %>% mutate(protein = prot_col)

# Resolve differing intensities
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

unify_intensity <- function(index, intensity_frame) {
  possible <- intensity_frame[index, ][!is.na(intensity_frame[index, ])]
  return(Mode(possible))
}

intensities <- lapply(file_names, function(x) {
  current <- joined %>% select(grep(x, colnames(.)))
  return(lapply(seq_along(1:dim(current)[1]), unify_intensity,
                intensity_frame = current) %>% unlist())
  }) %>% `names<-`(file_names)

final_frame <- tibble(protein = prot_col, ion = joined$ion) %>%
  bind_cols(as_tibble(intensities))

write_delim(final_frame, args$output, delim = "\t")
