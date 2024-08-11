PALETTE <- "PNWColors::Shuksan2"
library("glue")
library("seqinr")
library("gridExtra")
library("paletteer")
library("tidyverse")
library("reticulate")
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  tools <- "/home/shannc/Bio_SDD/tools"
} else {
  wd <- "/home/shannc/workflow"
  tools <- "/home/shannc/workflow/tools"
}

CATPUCCIN_LATTE <- list(
  rosewater = "#F2E9E1",
  flamingo = "#F2CDCD",
  pink = "#F4B8E4",
  mauve = "#CA9EE6",
  red = "#E78284",
  maroon = "#EA999C",
  peach = "#EF9F76",
  yellow = "#E5C890",
  green = "#A6D189",
  teal = "#81C8BE",
  sky = "#99D1DB",
  sapphire = "#85C1DC",
  blue = "#8CAAEE",
  lavender = "#BABBF1",
  text = "#4C4F69",
  subtext1 = "#5C5F77",
  subtext0 = "#6C6F85",
  overlay2 = "#7C7F93",
  overlay1 = "#8C8FA1",
  overlay0 = "#9CA0B0",
  surface2 = "#ACB0BE",
  surface1 = "#BCC0CC",
  surface0 = "#CCD0DA",
  base = "#EFF1F5",
  mantle = "#E6E9EF",
  crust = "#DCE0E8"
)


GRAPHS <- list()

run <- "ptools"
RESULTS <- glue("{wd}/results")
bm_data <- glue("{wd}/data/benchmark/ProteomeTools")
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo"),
  go_tm_dir = glue("{wd}/data/reference/.go_texts"),
  mapping = glue("{wd}/results/ptools/Databases/seq-header_mappings.tsv")
)
# CHOSEN_PASS <- "1-First_pass"
CHOSEN_PASS <- "2-Second_pass"
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))


srm_pool <- "Thermo_SRM_Pool_34"
proteotypic_pool <- "TUM_first_pool_87"
missing_pool <- "TUM_third_pool_4"
srm <- read_tsv(glue("{bm_data}/SRMAtlas Set_.tsv")) |> filter(`Pool name` == srm_pool)
proteotypic <- read_tsv(glue("{bm_data}/Proteotypic Set_.tsv")) |> filter(`Pool name` == proteotypic_pool)
missing <- read_tsv(glue("{bm_data}/MissingGeneSet.tsv")) |> filter(`Pool name` == missing_pool)



casanavo <- fasta2tb(glue("{wd}/results/ptools/Denovo/Casanovo/casanovo_normal.fasta")) |> distinct(sequence)
pepnet <- fasta2tb(glue("{wd}/results/ptools/Denovo/PepNet/pepnet_normal.fasta")) |> distinct(sequence)
smsnet <- fasta2tb(glue("{wd}/results/ptools/Denovo/SMSNet/smsnet_normal.fasta")) |> distinct(sequence)

all_denovo <- bind_rows(casanavo, pepnet, smsnet)
TRUE_PEP_TB <- bind_rows(missing, srm, proteotypic)


TRUE_PROT_TB <- fasta2tb(glue("{wd}/data/benchmark/ProteomeTools/TUM_third_pool_4_Thermo_SRM_Pool_34_TUM_first_pool_87.fasta")) |>
  mutate(header = map_chr(header, \(x) {
    str_remove(x, "^ACC:")
  }))


TRUE_PEPS <- TRUE_PEP_TB$Sequence |> map_chr(clean_peptide)
TRUE_PROT <- TRUE_PROT_TB$header

true_positive_denovo <- sum(TRUE_PEPS %in% all_denovo$sequence)


get_true_positives <- function(peptides) {
  tp <- sum(peptides %in% TRUE_PEPS)
  c(tp, tp / length(peptides))
}

reticulate::source_python(glue("{wd}/bin/helpers.py"))

main <- function(prefix, pass) {
  browser()
  result <- list(source = c(), n_true_positives = c(), true_positive_ratio = c())
  peptides <- read_tsv(glue("{RESULTS}/{prefix}/{pass}/percolator_peptide_map.tsv"))
  data <- read_tsv(glue("{RESULTS}/{prefix}/{pass}/{prefix}_all_wcoverage.tsv")) |>
    mutate(header = map_chr(header, \(x) {
      str_remove(x, "^ACC:")
    }))
  perc_proteins <- distinct(peptides, ProteinId, engine, .keep_all = TRUE)

  unique_peps <- flatten_by(data$unique_peptides, ";") |> unique()
  engines <- unique(peptides$engine)

  all_peps <- lapply(engines, \(x) {
    cur_filtered <- peptides |> filter(engine == x)
    cur_filtered$peptideIds %>% unique()
  }) %>%
    `names<-`(engines)
  # TODO: do this
  # all_headers <-
  all_peps$combined <- unique_peps
  match_comparison <- lapply(names(all_peps), \(x) {
    result$source <<- append(result$source, x)
    peps <- all_peps[[x]]
    headers <- NULL
    # headers <- all_headers[[x]]
    tp <- get_true_positives(peps, headers)
    result$n_true_positives <<- append(result$n_true_positives, tp[1])
    result$true_positive_ratio <<- append(result$true_positive_ratio, tp[2])
    find_matches(peps, TRUE_PEPS) |>
      as_tibble() |>
      mutate(source = x)
  }) |>
    bind_rows()
  list(metrics = as_tibble(result), tb = match_comparison)
}

outdir <- glue("{wd}/docs/figures/benchmark")
f1c <- glue("{outdir}/ptools_compare_1.tsv")
f1m <- glue("{outdir}/ptools_metrics_1.tsv")
f2c <- glue("{outdir}/ptools_compare_2.tsv")
f2m <- glue("{outdir}/ptools_metrics_2.tsv")

get_id_files <- function(comparison, metrics) {
  list(tb = read_tsv(comparison), metrics = read_tsv(metrics))
}


if (!file.exists(f1c)) {
  first <- main("ptools", "1-First_pass")
  write_tsv(first$metrics, f1m)
  write_tsv(first$tb, f1c)
} else {
  first <- get_id_files(f1c, f1m)
}

if (!file.exists(f2c)) {
  second <- main("ptools", "2-Second_pass")
  write_tsv(second$metrics, f2m)
  write_tsv(second$tb, f2c)
  second <- get_id_files(f2c, f2m)
}

main("ptools", "1-First_pass")
