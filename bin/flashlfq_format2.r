library(tidyverse)
library(optparse)

flashlfq_header <- c(
  "File Name", "Scan Retention Time", "Precursor Charge",
  "Base Sequence", "Full Sequence", "Peptide Monoisotopic Mass",
  "Protein Accession"
)
old_names <- c("file", "retensionTime", "precursorCharge",
               "base_peptide", "peptide", "mw", "protein")
parser <- OptionParser()
parser <- add_option(parser, c("-p, --path"),
  type = "character",
  help = "Path to formatted files"
)
parser <- add_option(parser, c("-o", "--output"),
  type = "character",
  help = "Output file name"
)
mapping <- read.delim(args$mapping, sep = "\t")
args <- parse_args(parser)
files <- paste0(args$path, list.files(args$path))

# for testing
## path <- "/home/sc31/peps/"
## files <- paste0(path, list.files(path))

all_engines <- lapply(files, function(x) {
  current <- read.delim(x, sep = "\t") %>%
    select(old_names) %>%
    rename_with(~flashlfq_header, all_of(old_names)) %>%
    as_tibble()
}) %>% bind_rows()

all_engines <- all_engines %>%
  filter(!(is.na(`Peptide Monoisotopic Mass`))) %>%
  mutate(`Base Sequence` = unlist(lapply(`Base Sequence`, gsub, pattern = "X", replacement = ""))) %>%
  distinct(`Base Sequence`, .keep_all = TRUE)

write_delim(all_engines, args$output, delim = "\t")
