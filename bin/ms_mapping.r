library(tidyverse)
library(MSnbase)
## Read a mzML file and generate a mapping of msms scans to their
##  precursor intensities, charges and retention time for use with
##  quantification

args <- commandArgs(trailingOnly = TRUE)

## args <- c("~/CiCs1.mzML", "~/CiCs_metrics.tsv") # Testing
input <- args[1]
output <- args[2]
run_name <- gsub("\\..*", "", input)

mzml <- header(openMSfile(input)) %>%
  as_tibble()
mzml <- mzml %>%
  mutate(scanNum = paste0(run_name, ".", acquisitionNum)) %>%
  select(c(scanNum, msLevel, retentionTime, precursorCharge,
           precursorIntensity, precursorMZ, totIonCurrent,
           peaksCount, basePeakMZ))

write_delim(mzml, output, delim = "\t")
