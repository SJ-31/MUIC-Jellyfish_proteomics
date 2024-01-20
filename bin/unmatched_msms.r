library(MSnbase)
library(tidyverse)
library(optparse)
library(glue)

# Take psm2combined pep file from all engines, concatenate them,
# Take all scans from all engines, concatenate to keep only scan and file, also
# remove duplicates

filter_msms <- function(file_name, msms_path, df, metric_df) {
  # Read in an mzML file of "file_name", filter it to leave ms 1 & 2 spectra
  # not found in "df", then write the filtered file
  df <- filter(df, file == file_name)
  n_matched <- df$scan %>%
    unique() %>%
    length()
  imported_ms <- readMSData(msms_path, mode = "onDisk")
  header <- header(imported_ms)
  n_all_msms <- (header %>% filter(msLevel == 2) %>% dim())[1]
  unmatched_array <- !(header$acquisitionNum %in% as.numeric(df$scan))
  unmatched_mzml <- imported_ms[unmatched_array]
  metrics <- tibble(
    file = file_name, n_msms = n_all_msms,
    n_matched_msms = n_matched,
    n_unmatched = n_all_msms - n_matched,
    percent_matched = n_matched / n_all_msms
  )
  writeMSData(unmatched_mzml, glue("filtered_{file_name}.mzML"))
  return(metrics)
}

filter_from_mzML <- function(msms_files, scan_df) {
  msms_files <- msms_files %>%
    `names<-`(gsub("\\..*", "", msms_files)) %>%
    as.list()
  metrics <- lapply(names(msms_files), function(x) {
    filter_msms(x, msms_files[[x]], scan_df)
  })
  return(metrics %>% bind_rows())
}


group_files <- function(file_list, type, selection, distinct_var) {
  # Read in either a psm2prot or mapped scan file
  if (type == "psm") {
    files <- lapply(file_list, function(x) {
      changed <- read_csv(x) %>%
        select(all_of(selection)) %>%
        filter(Is_decoy == FALSE) %>%
        mutate(peptide = unlist(lapply(peptide, cleanPeptide)))
      return(changed)
    })
  } else {
    files <- lapply(file_list, function(x) {
      return(read_tsv(x) %>% select(all_of(selection)))
    })
  }
  combined <- bind_rows(files)
  return(combined %>% distinct(!!as.symbol(distinct_var), .keep_all = TRUE))
}

join_psms_scans <- function(psm_list, scan_list, pep_thresh) {
  write(glue("PEP threshold: {pep_thresh}", stdout()))
  scans <- group_files(
    scan_list, "scan", c("file", "scan", "base_peptide"),
    "base_peptide"
  )
  psms <- group_files(psm_list, "psm", c("PEP", "peptide", "Is_decoy"), "peptide")
  join <- inner_join(scans, psms, by = join_by(x$base_peptide == y$peptide)) %>%
    filter(PEP <= pep_thresh) %>%
    mutate(scan = unlist(gsub(".*\\.", "", scan)))
  return(join)
}

main <- function(args) {
  ## Extract ms/ms spectra that haven't matched to anything
  psms <- glue("{args$psm_path}/{list.files(path = args$psm_path)}")
  scans <- glue("{args$scan_path}/{list.files(path = args$scan_path)}")
  mzmls <- list.files(args$mzML_path, pattern = "*.mzML")
  joined_psms <- join_psms_scans(
    psm_list = psms, scan_list = scans,
    pep_thresh = args$pep_thresh
  )
  metrics <- filter_from_mzML(mzmls, joined_psms)
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  parser <- OptionParser()
  parser <- add_option(parser, c("-m", "--psm_path"))
  parser <- add_option(parser, c("-s", "--scan_path"))
  parser <- add_option(parser, c("-p", "--pep_thresh"))
  parser <- add_option(parser, c("-z", "--mzML_path"))
  parser <- add_option(parser, c("-r", "--r_source"))
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  main(args)
}
