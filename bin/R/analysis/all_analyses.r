library("gt")
library("typed")
library("tidyverse")
library("ggplot2")
library("ggridges")
library("paletteer")
library("Peptides")
library("glue")

save <- function(to_save, outdir) {
  byType <- function(name, object) {
    if ("gg" %in% class(object) || "grob" %in% class(object)) {
      ggsave(glue("{outdir}/{name}.png"), object, width = 10, height = 10)
    } else if ("gt_tbl" %in% class(object)) {
      gtsave(object, glue("{outdir}/{name}.html"))
      gtsave(object, glue("{outdir}/{name}.tex"))
    } else if ("tbl_df" %in% class(object)) {
      write_tsv(object, glue("{outdir}/{name}.tsv"))
    }
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  lapply(names(to_save), \(x) {
    byType(x, to_save[[x]])
  })
}


M <- list()
if (str_detect(getwd(), "Bio_SDD")) {
  M$wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  M$env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  M$tools <- "/home/shannc/Bio_SDD/tools"
} else {
  M$wd <- "/home/shannc/workflow"
  M$env <- "/home/shannc/anaconda3/envs/reticulate"
  M$tools <- "/home/shannc/workflow/tools"
}
M$chosen_pass <- "1-First_pass"
M$path <- glue("{M$wd}/results/C_indra")
M$fdr <- 0.05
M$outdir <- glue("{M$path}/Analysis")
M$prottrans_embd <- glue("{M$outdir}/Embeddings_prottrans/embeddings.hdf5")
M$prottrans_dist <- glue("{M$outdir}/Embeddings_prottrans/distances.hdf5")
M$percolator_all_path <- glue("{M$path}/{M$chosen_pass}/percolator_all.tsv")
M$unmatched_path <- glue("{M$path}/{M$chosen_pass}/Quantify/Unmatched/unmatched_peptides.tsv")
M$seq_map_path <- glue("{M$path}/Databases/seq-header_mappings.tsv")
M$aligned_peptides_path <- glue("{M$path}/{M$chosen_pass}/aligned_peptides.tsv")
M$sample_name <- "C_indra"
M$r_source <- glue("{M$wd}/bin/R")
M$python_source <- glue("{M$wd}/bin")
M$embd_type <- "protein"
M$uniprot_data_dir <- glue("{M$wd}/data/protein_databases/comparison_taxa")
M$combined_results <- glue("{M$path}/{M$chosen_pass}/{M$sample_name}_all_wcoverage.tsv")
M$peptide_map_path <- glue("{M$path}/{M$chosen_pass}/percolator_peptide_map.tsv")
M$ontologizer_path <- glue("{M$outdir}/Ontologizer")
M$embedding_path <- glue("{M$wd}/data/reference/go_embedded.npz")
M$ontologizer_exec <- glue("{M$tools}/Ontologizer.jar")
M$orgdb_path <- glue("{M$wd}/tests/testthat/output/org.Cindrasaksajiae.eg.db")
M$go_path <- glue("{M$wd}/data/reference/go.obo")
M$go_slim_path <- glue("{M$wd}/data/reference/goslim_generic.obo")
M$go_tm_dir <- glue("{M$wd}/data/reference/.go_texts")
M$eggnog_cols <- c(
  "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
  "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
  "CAZy", "BiGG_Reaction", "PFAMs"
)
args <- list(
  r_source = M$r_source, python_source = M$python_source,
  go_info = glue("{M$outdir}/all_go_info.tsv"),
  go_tm_dir = glue("{M$wd}/data/reference/.go_texts")
)

source(glue("{M$r_source}/helpers.r"))
source(glue("{M$r_source}/GO_helpers.r"))
source(glue("{M$r_source}/cluster_helpers.r"))
source(glue("{M$r_source}/analysis/prepare_embeddings.r"))
source(glue("{M$r_source}/DR_helpers.r"))
source(glue("{M$r_source}/analysis/metric_functions.r"))

if (!dir.exists(M$outdir)) {
  dir.create(M$outdir)
} else if (!dir.exists(glue("{M$outdir}/Figures"))) {
  dir.create(glue("{M$outdir}/Figures"))
}

M$run <- get_run("C_indra", M$path)
M$taxa_tb <- read_tsv(glue("{M$path}/{M$chosen_pass}/{M$sample_name}_taxonomy.tsv"))
M$data <- read_tsv(glue("{M$outdir}/Aggregated/{M$sample_name}_all_wcategory.tsv"))
M$repr <- read_tsv(glue("{M$outdir}/{M$sample_name}_all_representatives.tsv"))
M$alignments <- get_alignment_data(M$path, M$chosen_pass)

GET_GO <- FALSE
if (GET_GO) {
  d <- get_go_data(M$combined_results,
    onto_path = M$ontologizer_path
  )
  if (!file.exists(glue("{M$outdir}/all_go_info.tsv"))) {
    get_go_vec(d$sample_tb, go_column = "GO_IDs", unique = TRUE) |>
      go_info_tb() |>
      write_tsv(glue("{M$outdir}/all_go_info.tsv"))
  }
  if (!file.exists(glue("{M$outdir}/all_go_slims.tsv"))) {
    get_go_vec(d$sample_tb, go_column = "GO_slims", unique = TRUE) |>
      go_info_tb() |>
      write_tsv(glue("{M$outdir}/all_go_slims.tsv"))
  }
}

sem_matrix_file <- glue("{M$outdir}/sem_matrices.hdf5")
if (!file.exists(sem_matrix_file)) {
  source(glue("{M$r_source}/analysis/save_sim.r"))
}

PLOT_GO <- FALSE
if (PLOT_GO) {
  go_map_ref <- read_tsv(glue("{M$wd}/data/reference/go_map_generic.tsv"))
  M$go_abbrev_file <- glue("{M$wd}/data/reference/go_abbrevs.json")
  source(glue("{M$r_source}/GO_text_mining_helpers.r"))
  clouds <- special_go_clouds(M$data)
  lapply(names(clouds), \(x) {
    ggsave(glue("{M$outdir}/Figures/{x}_go_word_cloud.png"), clouds[[x]], width = 15, height = 15)
  })
}

# NOTE: Name each element of the list by file it will be saved to
# ------------------------
# TODO
#' @details
#' Misc includes the following
#' - Trends in missing quantification
#' - Demo that grouping property is not conserved
#' - Isolate cases where peptides are longer than proteins they were matched to
#'   by blast
SOURCED <- TRUE
