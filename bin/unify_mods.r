library(tidyverse)
file <- commandArgs(trailingOnly = TRUE)
# Unify n-terminal modifications, remove indication of fixed modifications
nA <- "n[42.0106]"
mOx <- "M[15.9949]"

identipy_samples <- c("acetyl-ALPEHEK", "SSLoxMNSGR")
msfragger_samples <- c("R.M[15.9949]NM[15.9949]DNHRTNINDR.R", "-.SILAAEVYADLTFM[15.9949]AK.P",
    "K.TVYC[57.0215]EEKFNK.I")
comet_samples <- c("K.M[15.9949]RLM[15.9949]HR.R")

my_insert <- function(insert, pos, string) {
  # Insert a string into another string
  splits <- strsplit(string, "")[[1]]
  if (pos < 0) {
    pos <- length(splits) + 1 - abs(pos)
  }
  begin <- splits[1:pos-1]
  begin <- begin[!is.na(begin)]
  end <- splits[pos:length(splits)]
  end <- end[!is.na(end)]
  return(paste0(c(begin, insert, end), collapse = ""))
}

identipy_mods <- list("cam" = "C", "oxM" = mOx, "acetyl-" = nA)

change_mods <- function(modified_seq, mod_map) {
    for (i in seq_along(mod_map)) {
      modified_seq <- gsub(names(mod_map)[i], mod_map[i], modified_seq)
    }
    return(modified_seq)
}

identipy_mods(identipy_samples[2], identipy_mods)
