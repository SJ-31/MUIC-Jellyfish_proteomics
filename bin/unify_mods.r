library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
engine <- args[2]
file <- "./test_manifest_identipy.tsv"
engine <- "comet"
# Unify n-terminal modifications, remove indication of fixed modifications
nA <- "n[42.0106]"
mOx <- "M[15.9949]"
metamorpheus_mods <- list("[Common Fixed:Carbamidomethyl on C]" = "C",
                          "[Common Variable:Oxidation on M]" = mOx)
identipy_mods <- list("cam" = "C", "oxM" = "M", "acetyl-" = nA)
fragger_mods <- list("C[57.0215]" = "C")

change_mods <- function(modified_seq, mod_map) {
    for (i in seq_along(mod_map)) {
      modified_seq <- gsub(names(mod_map)[i], mod_map[i], modified_seq)
    }
    return(modified_seq)
}
ipy <- read.delim(file, sep = "\t") %>% as_tibble()
ipy <- ipy %>% mutate(Peptide = change_mods(Modified.sequence, identipy_mods))
ipy$Peptide

## if (engine == "comet") {

## } else if (engine == "msfragger") {

## } else if (engine == "maxquant") {

## } else if (engine == "identipy") {
##   file <- read.table(sep

## } else if (engine == "metamorpheus") {

## }
