library("gt")
library("plotly")
library("typed")
library("tidyverse")
library("ggplot2")
library("ggridges")
library("paletteer")
library("Peptides")
library("glue")


save <- function(to_save, outdir) {
  by_type <- function(name, object) {
    if ("gg" %in% class(object) || "grob" %in% class(object)) {
      width <- ifelse(is.null(attr(object, "width")), 10, attr(object, "width"))
      height <- ifelse(is.null(attr(object, "height")), 10, attr(object, "height"))
      ggsave(glue("{outdir}/{name}.png"), object, width = width, height = height)
    } else if ("gt_tbl" %in% class(object)) {
      gtsave(object, glue("{outdir}/{name}.html"))
      gtsave(object, glue("{outdir}/{name}.tex"))
    } else if ("tbl_df" %in% class(object)) {
      write_tsv(object, glue("{outdir}/{name}.tsv"))
    } else if ("plotly" %in% class(object) && "htmlwidget" %in% class(object)) {
      plotly::save_image(object, glue("{outdir}/{name}.svg"), width = 1000, height = 800)
    } else if (object == 0) {
      cat("", file = glue("{outdir}/{name}.txt"))
    }
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  lapply(names(to_save), \(x) {
    by_type(x, to_save[[x]])
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

# ----------------------------------------
#' Main entry point to choose paths
M$prefixes <- list("C_indra", "C_indra.msconvert", "C_indra.calibrated", "ND_C_indra")
M$chosen_prefix <- prefixes[[1]]
M$chosen_pass <- "2-Second_pass"
M$chosen_path <- glue("{M$wd}/results/C_indra")


M$path <- glue("{M$wd}/results/C_indra")
M$cpath <- glue("{M$wd}/results/C_indra.calibrated")
M$mpath <- glue("{M$wd}/results/C_indra.msconvert")
M$ndpath <- glue("{M$wd}/results/ND_C_indra")
M$all_paths <- list(M$path, M$mpath, M$cpath, M$ndpath)
M$fdr <- 0.05
M$outdir <- glue("{M$wd}/docs/figures")
M$prottrans_embd <- glue("{M$outdir}/Embeddings_prottrans/embeddings.hdf5")
M$prottrans_dist <- glue("{M$outdir}/Embeddings_prottrans/distances.hdf5")
M$percolator_all_path <- glue("{M$chosen_path}/{M$chosen_pass}/percolator_all.tsv")
M$seq_map_path <- glue("{M$chosen_path}/Databases/seq-header_mappings.tsv")
M$aligned_peptides_path <- glue("{M$chosen_path}/{M$chosen_pass}/aligned_peptides.tsv")
M$r_source <- glue("{M$wd}/bin/R")
M$python_source <- glue("{M$wd}/bin")
M$embd_type <- "protein"
M$uniprot_data_dir <- glue("{M$wd}/data/protein_databases/comparison_taxa")
M$data_path <- glue("{M$chosen_path}/{M$chosen_pass}/{M$chosen_prefix}_all_wcoverage.tsv")
M$peptide_map_path <- glue("{M$chosen_path}/{M$chosen_pass}/percolator_peptide_map.tsv")
M$ontologizer_path <- glue("{M$chosen_path}/Analysis/Ontologizer")
M$embedding_path <- glue("{M$wd}/data/reference/go_embedded.npz")
M$ontologizer_exec <- glue("{M$tools}/Ontologizer.jar")
M$orgdb_path <- glue("{M$chosen_path}/Analysis/{M$chosen_prefix}_org.db")
M$lfq <- read_tsv(glue("{M$chosen_path}/{M$chosen_pass}/lfq_all.tsv"))
M$go_path <- glue("{M$wd}/data/reference/go.obo")
M$go_slim_path <- glue("{M$wd}/data/reference/goslim_generic.obo")
M$go_tm_dir <- glue("{M$wd}/data/reference/.go_texts")
M$aq_reformat_path <- glue("{M$chosen_path}/{M$chosen_pass}/Quantify/directlfq.aq_reformat.tsv")
M$eggnog_cols <- c(
  "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
  "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
  "CAZy", "BiGG_Reaction", "PFAMs"
)
M$mapped_scan_path <- glue("{M$chosen_path}/{M$chosen_pass}/all_mapped_scans.tsv")
if (!file.exists(M$mapped_scan_path)) {
  tb <- lapply(c("Quantify", "Open_search"), \(x) {
    list.files(glue("{M$chosen_path}/{M$chosen_pass}/{x}/Mapped_scans"),
      pattern = "*.tsv",
      full.names = TRUE
    )
  }) |>
    unlist() |>
    lapply(read_tsv) |>
    bind_rows()
  write_tsv(tb, M$mapped_scan_path)
}

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
}

# Note that you must match on headers
M$run <- get_run(M$chosen_prefix, M$chosen_path)
M$taxa_tb <- read_tsv(glue("{M$chosen_path}/{M$chosen_pass}/{M$chosen_prefix}_taxonomy.tsv"))
if (str_detect(M$chosen_pass, "Second")) {
  M$data <- M$run$second
} else {
  M$data <- M$run$first
}
if (file.exists(glue("{M$chosen_path}/Analysis/{M$chosen_prefix}_all_wcog.tsv"))) {
  M$data_w_cat <- read_tsv(glue("{M$chosen_path}/Analysis/{M$chosen_prefix}_all_wcog.tsv"))
}
M$alignments <- get_alignment_data(M$chosen_path, M$chosen_pass)
M$taxa_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

GET_GO <- FALSE
if (GET_GO) {
  d <- get_go_data(M$data_path,
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

PLOT_GO <- FALSE
if (PLOT_GO) {
  go_map_ref <- read_tsv(glue("{M$wd}/data/reference/go_map_generic.tsv"))
  source(glue("{M$r_source}/GO_text_mining_helpers.r"))
  clouds <- special_go_clouds(M$data)
  lapply(names(clouds), \(x) {
    ggsave(glue("{M$outdir}/Figures/{x}_go_word_cloud.png"), clouds[[x]], width = 15, height = 15)
  })
}

# NOTE: Name each element of the list by file it will be saved to
# ------------------------
#' @details
#' Misc includes the following
#' - Trends in missing quantification
#' - Demo that grouping property is not conserved
#' - Isolate cases where peptides are longer than proteins they were matched to
#'   by blast
SOURCED <- TRUE
