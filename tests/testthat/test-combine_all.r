library(glue)
library(testthat)
library(ggplot2)
source("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/R/combine_all.r")

setwd("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/debug/combine")

args <- list(eggnog = "C_indra_eggnog_matched.tsv",
             interpro = "C_indra_interpro_matched.tsv",
             downloads = "C_indra_downloads_anno-3.tsv",
             sort_mods = TRUE,
             empai = TRUE,
             pfam2go = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/pfam2go",
             interpro2go = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/interpro2go",
             pfam_db = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Databases/pfam_entries.tsv",
             is_denovo = "false",
             fdr = 0.05,
             pep_thresh = 1,
             output = "C_indra_all.tsv",
             directlfq = "sorted_directlfq.tsv",
             maxlfq = "max_lfq.tsv",
             flashlfq = "sorted_flashlfq.tsv",
             r_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin/R", python_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/GO_helpers.r"))
results <- main(args)
all <- results$all
found <- results$f

getDuplicates <- function(tb) {
  dupes <- tb %>% group_by(ProteinId) %>% filter(n() > 1)
  print(glue("N duplicates: {nrow(dupes)}"))
  return(dupes)
}


checkNa <- function(matrx, name) {
  print(glue("testing {name}"))
  expect_equal(matrx[["Anno_method"]], 0)
  expect_equal(matrx[["ID_method"]], 0)
  expect_equal(matrx[["header"]], 0)
  expect_equal(matrx[["ProteinId"]], 0)
  expect_equal(matrx[["q.value"]], 0)
  print(glue("testing {name} DONE"))
}

strangeNas <- function(args) {
  eggnog <- loadAndMerge(args$eggnog_anno, args$eggnog_meta)
  interpro <- loadAndMerge(args$interpro_anno, args$interpro_meta)
  downloads <- loadAndMerge(args$download_anno, args$download_meta)
  eggnog_nas <- colSums(is.na(eggnog))
  checkNa(eggnog_nas, "eggnog")
  interpro_nas <- colSums(is.na(interpro))
  checkNa(interpro_nas, "interpro")
  downloads_nas <- colSums(is.na(downloads))
  checkNa(downloads_nas, "downloads")
  browser()
  print("temp")
}

testCoverage <- function(df) {
  source(glue("{args$r_source}/protein_coverage.r"))
  start <- Sys.time()
  covered <- coverageCalc(df)
  end <- Sys.time()
  print(glue("Coverage took {end-start}"))
  return(covered)
}


testResolveAlignment <- function(wf) {
  source(glue("{args$r_source}/protein_coverage.r"))
  expect_equal(resolveAlignment(
    "--AGFS--", "DSAGFS--"
  ), "DSAGFS--")
  expect_equal(resolveAlignment(
    "--[AG]FS--", "--GFS--"
  ), "--[AG]FS--")
  expect_equal(resolveAlignment(
    "--[AG]FS--", "--DFS--"
  ), "--[AGD]FS--")
  expect_equal(resolveAlignment(
    "--AGFS--", "--AG----"
  ), "--AGFS--")
  expect_equal(resolveAlignment(
    "--AGFS--", "--DF----"
  ), "--[AD][GF]FS--")
}


testResolveAlignment2 <- function(wf) {
  source(glue("{args$r_source}/protein_coverage.r"))
  expect_equal(resolveAlignment2(
    "--AGFS--", "DSAGFS--"
  ), list("DSAGFS--"))
  expect_equal(resolveAlignment2(
    "--AGFS--", "--AFSD--"
  ), list("--AGFS--", "--AFSD--"))
  expect_equal(resolveAlignment2(
    list("--AGFS--", "--AFSD--"), "--AGFSFF"
  ), list("--AGFS--", "--AFSD--", "--AFSDFF"))
}

## testResolveAlignment()

sample <- dplyr::slice(all, 1:10)
s <- testCoverage(sample)

WANTED_TYPES <- c("list", "character", "logical", "numeric")
UNWANTED_TYPES <- c("tbl", "data.frame", "tbl_df", "function")

checkWanted <- function(x) {
  evaluated <- get(x)
  what_is <- class(evaluated)
  if (any(what_is %in% UNWANTED_TYPES)) {
    return(FALSE)
  }
  if (purrr::pluck_depth(evaluated) > 1) {
    any(lapply(seq_along(evaluated), \(x) {
      any(evaluated[[x]] %in%
            UNWANTED_TYPES)
    }) %>% unlist())
    return(FALSE)
  }
  return(any(what_is %in% WANTED_TYPES))
}

kept_vars <- ls() %>% keep(., checkWanted)
with_vals <- lapply(kept_vars, get) %>% `names<-`(kept_vars)
