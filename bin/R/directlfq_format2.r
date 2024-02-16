library(tidyverse)
library(optparse)
library(glue)


distinctProt <- function(protein_col) {
  unlist(lapply(protein_col, \(x) {
    unique(str_split_1(x, ";")) %>% paste0(collapse = ";")
  }))
}


medIntensity <- function(intensity_col) {
  unlist(lapply(intensity_col, \(x) {
    unique(str_split_1(x, ";")) %>%
      as.double() %>%
      median()
  }))
}


file_pivot <- function(psm_df) {
  sort_ions <- psm_df %>%
    filter(!is.na(protein)) %>%
    group_by(ion, file) %>%
    mutate_at(., vars(-group_cols()), paste0, collapse = ";") %>%
    distinct() %>%
    ungroup()
  sort_ions <- sort_ions %>% mutate(
    protein = distinctProt(protein),
    precursorIntensity = medIntensity(precursorIntensity)
  )
  pivoted <- sort_ions %>% pivot_wider(
    names_from = file,
    values_from = precursorIntensity
  )
  return(pivoted)
}

main <- function(args) {
  mapping <- read.delim(args$mapping, sep = "\t")
  file_names <- mapping$scanNum %>%
    lapply(., gsub, pattern = "\\..*", replacement = "") %>%
    unlist(use.names = FALSE) %>%
    unique()
  rm(mapping)

  files <- paste0(args$path, "/", list.files(args$path))
  all_engines <- lapply(files, function(x) {
    df <- read.delim(x, sep = "\t") %>%
      as_tibble() %>%
      rename(ion = base_peptide) %>%
      select(c("file", "protein", "ion", "precursorIntensity")) %>%
      file_pivot()
    return(df)
  })

  # Merge matches between engines
  joined <- Reduce(function(x, y) {
    full_join(x, y, by = "ion")
  }, all_engines)


  # Group up proteins that share an ion
  prot <- joined %>% select(grep("protein", colnames(joined)))
  prot_col <- lapply(seq_along(1:dim(prot)[1]), function(x) {
    filtered <- prot[x, ][grepl("[a-zA-Z1-9]+", prot[x, ])]
    return(paste0(filtered, collapse = ";"))
  }) %>% unlist()
  joined <- joined %>% mutate(protein = prot_col)

  intensities <- lapply(file_names, \(x) {
    current <- select(joined, grep(x, colnames(joined)))
    medians <- apply(current, 1, median, na.rm = TRUE)
  }) %>% `names<-`(file_names)

  final_frame <- tibble(protein = prot_col, ion = joined$ion) %>%
    bind_cols(as_tibble(intensities)) %>%
    filter(protein != "") %>%
    mutate_all(~ replace(., is.na(.), 0))

  final_frame <- final_frame %>%
    group_by(ion) %>%
    mutate_at(., vars(-group_cols()), paste0, collapse = ";") %>%
    distinct() %>%
    ungroup() %>%
    mutate(protein = distinctProt(protein)) %>%
    mutate(across(-c("protein", "ion"), medIntensity))
  return(final_frame)
}


if (sys.nframe() == 0) { # Won't run if the script is being sourced
  parser <- OptionParser()
  parser <- add_option(parser, c("-p", "--path"),
    type = "character",
    help = "Path to formatted files"
  )
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  parser <- add_option(parser, c("-m", "--mapping"),
    type = "character",
    help = "Mapping file name"
  )
  args <- parse_args(parser)
  f <- main(args)
  write_delim(f, args$output, delim = "\t")
}
