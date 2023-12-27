library(testthat)
library(glue)
bin <- "../../bin"
source(glue("{bin}/get_scan_num.r"))

file_list <- list( # For tests
  comet = "../../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv",
  identipy = "../../results/test_manifest/1-First_pass/Percolator/identipy_percolator_psms.tsv", maxquant = "../../results/test_manifest/1-First_pass/MaxQuant/maxquant_all_pins.temp", msfragger = "../../results/test_manifest/1-First_pass/Percolator/msfragger_percolator_psms.tsv",
  metamorpheus = "../../results/test_manifest/1-First_pass/Metamorpheus/metamorpheus_AllPSMs.psmtsv", tide = "../../results/test_manifest/1-First_pass/Tide/tide-search.target.txt"
)
chosen <- "metamorpheus" # for testing
mapping <- read.delim("../toy_results/CiCs_metrics.tsv", sep = "\t")
final <- read_engine_psms(file_list[[chosen]], chosen, mapping)
final
