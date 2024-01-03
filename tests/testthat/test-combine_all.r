source("./bin/combine_all.r")
dir <- "./results/jellyfish/1-First_pass"
args <- list(
  eggnog_anno = glue("{dir}/Unmatched/eggNOG/eggnog_anno-nO_D.tsv"),
  eggnog_meta = glue("{dir}/Unmatched/eggNOG/eggnog_meta-nO_D.tsv"),
  interpro_anno = glue("{dir}/Unmatched/InterPro/interpro_anno-nO_D.tsv"),
  interpro_meta = glue("{dir}/Unmatched/InterPro/interpro_meta-nO_D.tsv"),
  download_anno = glue("{dir}/Unmatched/Database-annotated/download_anno-nO_D.tsv"),
  download_meta = glue("{dir}/Unmatched/Database-annotated/download_meta-nO_D.tsv"),
  coverage = FALSE,
  sort_mods = FALSE,
  empai = FALSE,
  is_denovo = "true",
  directlfq = glue("{dir}/Quantify/sorted_directlfq.tsv"),
  flashlfq = glue("{dir}/Quantify/sorted_flashlfq.tsv"),
  output_anno = "./tests/testthat/output/combined-anno.tsv",
  output_meta = "./tests/testthat/output/combined-anno.tsv",
  r_source = "./bin/",
  fdr = 0.05,
  pep_thresh = 1
)

results <- main(args)
all <- results$all
write_lines(results$anno$ProteinId, "./tests/testthat/output/all_proteinids.txt")
found <- results$f
