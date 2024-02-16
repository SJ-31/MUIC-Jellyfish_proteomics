library(tidyverse)
input <- "./msmsScans.txt"

decoy_label <- function(maxquant_proteins) {
  split <- unlist(strsplit(maxquant_proteins, ";", fixed = TRUE))
  if (substr(split[1], 1, 4) == "REV_") {
    return(-1)
  }
  return(1)
}
new_names <- c(Peptide = "Sequence", ScanNr = "Scan.number")

msms <- read.delim(input, sep = "\t")
wanted <- c("SpecId", "Label", "ScanNr",

  "Peptide", "Proteins")
msms %>% colnames

# msms <- msms %>%
#   mutate(SpecId = paste0(Raw.file, ".", Scan.number)) %>%
#   mutate(Label = unlist(lapply(Proteins, decoy_label))) %>%
#   mutate(Proteins = unlist(lapply(Proteins, gsub,
#     pattern = ";",
#     replacement = "\t"
#   ))) %>%
