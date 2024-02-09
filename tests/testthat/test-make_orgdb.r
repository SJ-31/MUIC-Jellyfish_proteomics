library(glue)
source("./bin/make_orgdb.r")

args <- list(
  combined_annotations = "./results/jellyfish/1-First_pass/jellyfish_all.tsv",
  path = "./data/protein_databases/comparison_taxa",
  maintainer = "Shann Chongwattananukul <shann.cho@student.mahidol.edu>",
  author = "Shann Chongwattananukul <shann.cho@student.mahidol.edu>",
  tax_id = "1758319",
  genus = "Chironex",
  species = "indrasaksajiae"
)

dfs <- prepDf(args$combined_annotations, args$path)
## gos <- as_tibble(dfs$go)
## output <- "./tests/testthat/output/"
## unlink(glue("{output}/org.Cindrasaksajiae.eg.db"), recursive = TRUE, force = TRUE)
## makeOrgPackage(
##   go = dfs$go, info = dfs$info,
##   version = "0.1",
##   maintainer = args$author,
##   outputDir = output,
##   author = args$maintainer,
##   tax_id = args$tax_id,
##   genus = args$genus,
##   species = args$species,
##   goTable = "go",
##   verbose = TRUE
## )
