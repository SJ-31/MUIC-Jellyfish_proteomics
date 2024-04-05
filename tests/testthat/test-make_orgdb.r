library(glue)
source("./bin/R/make_orgdb.r")

args <- list(
  combined_annotations = "./results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv",
  maintainer = "Shann Chongwattananukul <shann.cho@student.mahidol.edu>",
  author = "Shann Chongwattananukul <shann.cho@student.mahidol.edu>",
  path_to_extra = "NONE",
  tax_id = "1758319",
  genus = "Chironex",
  species = "indrasaksajiae"
)

dfs <- prepDf(args$combined_annotations, args$path)
gos <- as_tibble(dfs$go)
output <- "./tests/testthat/output/"
unlink(glue("{output}/org.Cindrasaksajiae.eg.db"), recursive = TRUE, force = TRUE)
makeOrgPackage(
  go = dfs$go, info = dfs$info,
  version = "0.1",
  maintainer = args$author,
  outputDir = output,
  author = args$maintainer,
  tax_id = args$tax_id,
  genus = args$genus,
  species = args$species,
  goTable = "go",
  verbose = TRUE
)
