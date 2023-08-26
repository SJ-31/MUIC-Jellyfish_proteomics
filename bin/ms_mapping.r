library(tidyverse)
library(MSnbase)
## Read a mzML file and generate a mapping of msms scans to their
##  precursor intensities, charges and retension time for use with
##  quantification

args <- commandArgs(trailingOnly = TRUE)

## args <- c("../results/CiCs1.mzML", "../results/CiCs_metrics.tsv") # Testing
input <- args[1]
output <- args[2]

mzml <- readMSData(input)
spectra_list <- spectra(mzml)
run_name <- str_match(input, "(.*)\\..*")
format_scan <- function(spectra_name) {
  return(str_match(
    spectra_name,
    ".*\\.[A-Za-z]*0*([1-9]*.+)$"
  )[, 2])
}

scanIds <- names(spectra_list) %>%
  lapply(., format_scan) %>%
  paste0(run_name, ".", .)


mapped_scans <- data.frame(
  scanNum = scanIds,
  retensionTime = unlist(lapply(spectra_list, rtime), use.names = FALSE),
  precursorCharge = unlist(lapply(spectra_list, precursorCharge),
    use.names = FALSE
  ),
  precursorIntensity = unlist(lapply(spectra_list, precursorIntensity),
    use.names = FALSE
  )
)

write_delim(mapped_scans, output, delim = "\t")
