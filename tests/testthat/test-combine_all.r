library(testthat)
source("./bin/combine_all.r")
dir <- "./results/jellyfish/1-First_pass"

args <- list(
  eggnog = glue("{dir}/Unmatched/eggNOG/jelly_eggnog_matched.tsv"),
  interpro = glue("{dir}/Unmatched/InterPro/jelly_interpro_matched.tsv"),
  downloads =
    glue("{dir}/Unmatched/Database-annotated/jelly_downloads_anno-3.tsv"),
  coverage = FALSE,
  sort_mods = FALSE,
  empai = FALSE,
  is_denovo = "true",
  directlfq = glue("{dir}/Quantify/sorted_directlfq.tsv"),
  flashlfq = glue("{dir}/Quantify/sorted_flashlfq.tsv"),
  output = "./tests/testthat/output/combined-anno.tsv",
  r_source = "./bin/",
  fdr = 0.05,
  pep_thresh = 1
)

results <- main(args)
all <- results$all
write_tsv(all, "./tests/results/Combined/all_test.tsv")
write_lines(results$anno$ProteinId, "./tests/testthat/output/all_proteinids.txt")
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


strangeNas(args)
