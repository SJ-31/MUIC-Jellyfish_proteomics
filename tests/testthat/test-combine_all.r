library(glue)
library(testthat)
library(ggplot2)
source("./bin/combine_all.r")
dir_name <- "./results/jellyfish/1-First_pass"
tests <- "./tests/results/"

args <- list(
  eggnog = glue::glue("{dir_name}/Unmatched/eggNOG/jellyfish_eggnog_matched.tsv"),
  interpro = glue::glue("{dir_name}/Unmatched/InterPro/jellyfish_interpro_matched.tsv"),
  downloads =
    glue::glue("{dir_name}/Unmatched/Database-annotated/jellyfish_downloads_anno-3.tsv"),
  coverage = FALSE,
  sort_mods = TRUE,
  empai = TRUE,
  is_denovo = "true",
  pfam2go = "./results/jellyfish/Databases/pfam2go",
  interpro2go = "./results/jellyfish/Databases/interpro2go",
  pfam_db = "./results/jellyfish/Databases/pfam_entries.tsv",
  directlfq = glue::glue("{dir_name}/Quantify/sorted_directlfq.tsv"),
  flashlfq = glue::glue("{dir_name}/Quantify/sorted_flashlfq.tsv"),
  output = "./tests/testthat/output/combined-anno.tsv",
  r_source = "./bin/",
  fdr = 0.05,
  pep_thresh = 1
)
results <- main(args)
all <- results$all
## write_tsv(all, "./tests/results/Combined/all_test.tsv")
## write_lines(results$anno$ProteinId, "./tests/testthat/output/all_proteinids.txt")
found <- results$f


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

## test <- s[1, ]
## algn <- test[["alignment"]] %>% splitWithGroups()
## seq <- test[["seq"]] %>% str_split_1("")
## mid <- c()
## seq_new <- c()
## algn_new <- c()
## lapply(seq_along(seq), \(x) {
##   if (algn[x] == seq[x]) {
##     mid <<- c(mid, "|")
##     algn_new <<- c(algn_new, algn[x])
##     seq_new <<- c(seq_new, seq[x])
##   } else if (grepl("\\[", algn[x])) {
##     residues <- str_extract(algn[x], "\\[(.*)\\]", group = 1)
##     midpoint <- median(seq_len(residues))
##     if (str_length(algn[x]) %% 2 == 0) {
##       midpoint <- ceiling(midpoint)
##       residues <- R.utils::insert(residues, midpoint, "|")
##     }
##     padding <- str_dup(" ", (length(residues) - 1) / 2)
##     algn_new <<- str_flatten(c("[", residues, "]"))
##     seq_new <<- c(seq_new, glue("{padding}{seq[x]}{padding}"))
##     mid <<- c(mid, glue("{padding}|{padding}"))
##   } else {
##     mid <<- c(mid, " ")
##     algn_new <<- c(algn_new, algn[x])
##     seq_new <<- c(seq_new, seq[x])
##   }
## })
## full_algn <- lapply(list(seq_new, mid, algn_new), str_flatten) %>% unlist()
