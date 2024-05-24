library(gt)
library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}
GET_GO <- TRUE
PLOT_GO <- FALSE
CHOSEN_PASS <- "1-First_pass"
PATH <- glue("{wd}/results/C_indra")
FDR <- 0.05
OUTDIR <- glue("{PATH}/Analysis")
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR)
} else if (!dir.exists(glue("{OUTDIR}/figures"))) {
  dir.create(glue("{OUTDIR}/figures"))
}
PROTTRANS_EMBD <- glue("{OUTDIR}/Embeddings_prottrans/embeddings.hdf5")
PROTTRANS_DIST <- glue("{OUTDIR}/Embeddings_prottrans/distances.hdf5")
saved <- glue("{PATH}")

# TODO: Change any directories that use `test` into valid ones
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  embd_type = "protein",
  sample_name = "C_indra",
  uniprot_data_dir = glue("{wd}/data/protein_databases/comparison_taxa"),
  combined_results = glue("{PATH}/{CHOSEN_PASS}/{args$sample_name}_all_wcoverage.tsv"),
  ontologizer_path = glue("{wd}/tests/nf-test-out/ontologizer/"),
  embedding_path = glue("{wd}/data/reference/go_embedded.npz"),
  orgdb_path = glue("{wd}/tests/testthat/output/org.Cindrasaksajiae.eg.db"),
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo"),
  go_tm_dir = glue("{wd}/data/reference/.go_texts")
)
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))

EGGNOG_COLS <- c(
  "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
  "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
  "CAZy", "BiGG_Reaction", "PFAMs"
)

run <- runData("C_indra", PATH)
if (GET_GO) {
  d <- goData(args$combined_results,
    onto_path = args$ontologizer_path
  )
}
if (PLOT_GO) {
  source(glue("{args$r_source}/GO_text_mining_helpers.r"), local = th)
  clouds <- specialGoClouds(d$sample_tb)
  lapply(names(clouds), \(x) {
    ggsave(glue("{OUTDIR}/figures/{x}_go_word_cloud.png"), clouds[[x]])
  })
}

# NOTE: Name each element of the list by file it will be saved to
GRAPHS <- list()
TABLES <- list()
SAVE <- FALSE
# ------------------------
# source(glue("{args$r_source}/analysis/subanalyses/general_metrics.r"))
# source(glue("{args$r_source}/analysis/subanalyses/PTM_analyses.r")) # TODO:
# source(glue("{args$r_source}/analysis/subanalyses/general_metrics.r"))

# source(glue("{args$r_source}/analysis/fgsea.r")) # DONE
# source(glue("{args$r_source}/analysis/subanalyses/comparison_with_previous.r")) # DONE
# source(glue("{args$r_source}/analysis/subanalyses/engine_category_bias.r")) # DONE
# source(glue("{args$r_source}/analysis/subanalyses/pass_differences.r")) # DONE




if (SAVE) {
  lapply(names(GRAPHS), \(x) {
    ggsave(glue("{OUTDIR}/figures/{x}"), GRAPHS[[x]])
  })
  lapply(names(TABLES), \(x) {
    gtsave(glue("{OUTDIR}/{x}"), TABLES[[x]])
  })
}


# --------------------------------------------------------
# Lineage metrics


# --------------------------------------------------------

# Grouping proteins by their participation in KEGG pathways
current <- run$first
pwy <- groupPathways(current)

# --------------------------------------------------------
# Investigating trends in missing quantification
current <- run$first
cq <- c("directlfq", "maxlfq", "flashlfq")
missing_quant_tests <- list()
for (q in cq) {
  noq <- current %>% filter(is.na(!!as.symbol(glue("{cq}_mean"))))
  hasq <- current %>% filter(!ProteinId %in% noq$ProteinId)
  missing_quant_tests[[q]] <- wilcox.test(noq$num_peps, hasq$num_peps,
    alternative = "l"
  )
}
capture.output(missing_quant_tests, file = glue("{OUTDIR}/missing_quantification_tests.txt"))
rm(noq)
rm(hasq)
# --------------------------------------------------------
