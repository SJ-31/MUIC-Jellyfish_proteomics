library(tidyverse)
library(optparse)
library(glue)

file_pivot <- function(psm_df) {
  return(psm_df %>%
    distinct(ion, .keep_all = TRUE) %>%
    pivot_wider(names_from = file, values_from = precursorIntensity))
}

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
mapping <- read.delim(args$mapping, sep = "\t")

file_names <- mapping$scanNum %>%
  lapply(., gsub, pattern = "\\..*", replacement = "") %>%
  unlist(use.names = FALSE) %>%
  unique()
rm(mapping)

files <- paste0(args$path, list.files(args$path))

all_engines <- lapply(files, function(x) {
  df <- read.delim(x, sep = "\t")  %>%
    as_tibble() %>%
    rename(ion = base_peptide) %>%
    select(c("file", "protein", "ion", "precursorIntensity")) %>%
    file_pivot()
  return(df)
})

# Merge matches between engines
joined <- Reduce(function(x, y) { full_join(x, y, by = "ion") }, all_engines)
rm(all_engines)

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
