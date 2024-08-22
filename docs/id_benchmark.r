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
TABLES <- list()

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
source(glue("{args$r_source}/analysis/all_analyses.r"))
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
  c(tp, tp / length(peptides), tp / length(TRUE_PEPS))
}

get_true_positives_prot <- function(proteins) {
  tp <- sum(proteins %in% TRUE_PROT)
  c(tp, tp / length(proteins), tp / length(TRUE_PROT))
}

reticulate::source_python(glue("{wd}/bin/helpers.py"))

main <- function(prefix, pass) {
  result <- list(
    source = c(), n_true_positives = c(),
    true_positive_ratio = c(),
    ratio_total = c()
  )
  result_prot <- list(
    source = c(), n_true_positives = c(),
    true_positive_ratio = c(),
    ratio_total = c()
  )
  peptides <- read_tsv(glue("{RESULTS}/{prefix}/{pass}/percolator_peptide_map.tsv"))

  data <- read_tsv(glue("{RESULTS}/{prefix}/{pass}/{prefix}_all_wcoverage.tsv")) |>
    mutate(header = map_chr(header, \(x) {
      str_remove(x, "^ACC:")
    }))
  seq_map <- read_tsv(glue("{RESULTS}/{prefix}/{pass}/seq-header_map_found.tsv")) |>
    mutate(header = map_chr(header, \(x) {
      str_remove(x, "^ACC:")
    }))
  perc_proteins <- distinct(peptides, ProteinId, engine, .keep_all = TRUE)

  obs_proteins <- perc_proteins |>
    filter(!grepl("D", ProteinId)) |>
    inner_join(seq_map, by = join_by(x$ProteinId == y$id))

  unique_peps <- flatten_by(data$unique_peptides, ";") |> unique()
  engines <- unique(peptides$engine)

  headers <- data$header |> unique()

  result_prot$source <- "combined"
  tp_prot <- sum(headers %in% TRUE_PROT)
  result_prot$n_true_positives <- tp_prot
  result_prot$ratio_total <- tp_prot / length(TRUE_PROT)
  result_prot$true_positive_ratio <- tp_prot / length(headers)

  all_peps <- lapply(engines, \(x) {
    cur_filtered <- peptides |> filter(engine == x)
    cur_filtered$peptideIds %>% unique()
  }) %>%
    `names<-`(engines)

  all_prot <- lapply(engines, \(x) {
    cur_filtered <- obs_proteins |> filter(engine == x)
    cur_filtered$header %>% unique()
  }) %>%
    `names<-`(engines)

  all_peps$combined <- unique_peps
  t_query <- lapply(names(all_peps), \(x) {
    result$source <<- append(result$source, x)
    peps <- all_peps[[x]]
    prots <- all_prot[[x]]

    tp <- get_true_positives(peps)
    result$n_true_positives <<- append(result$n_true_positives, tp[1])
    result$true_positive_ratio <<- append(result$true_positive_ratio, tp[2])
    result$ratio_total <<- append(result$ratio_total, tp[3])

    tp2 <- get_true_positives_prot(prots)
    result_prot$source <<- append(result$source, x)
    result_prot$n_true_positives <<- append(result_prot$n_true_positives, tp2[1])
    result_prot$true_positive_ratio <<- append(result_prot$true_positive_ratio, tp2[2])
    result_prot$ratio_total <<- append(result_prot$ratio_total, tp2[3])

    find_matches(TRUE_PEPS, peps) |> # Check if all true peptides are represented
      # by found peptides
      as_tibble() |>
      mutate(source = x)
  }) |>
    bind_rows()
  obs_query <- lapply(names(all_peps), \(x) {
    # Check if the found peptides deviate strongly from the true peptides
    find_matches(all_peps[[x]], TRUE_PEPS) |>
      as_tibble() |>
      mutate(source = x)
  }) |>
    bind_rows()
  list(
    metrics = as_tibble(result),
    true_as_query = t_query, obs_as_query = obs_query, pmetrics = as_tibble(result_prot)
  )
}

outdir <- glue("{wd}/docs/figures/benchmark")
f1tq <- glue("{outdir}/ptools_t_query_1.tsv")
f1oq <- glue("{outdir}/ptools_obs_query_1.tsv")
f1m <- glue("{outdir}/ptools_metrics_1.tsv")
f1pm <- glue("{outdir}/ptools_pmetrics_1.tsv")

f2tq <- glue("{outdir}/ptools_t_query_2.tsv")
f2oq <- glue("{outdir}/ptools_obs_query_2.tsv")
f2m <- glue("{outdir}/ptools_metrics_2.tsv")
f2pm <- glue("{outdir}/ptools_pmetrics_2.tsv")

get_id_files <- function(t_query, obs_query, metrics, pmetrics) {
  list(
    true_as_query = read_tsv(t_query), metrics = read_tsv(metrics),
    obs_as_query = read_tsv(obs_query), pmetrics = read_tsv(pmetrics)
  )
}


if (!file.exists(f1tq)) {
  first <- main("ptools", "1-First_pass")
  write_tsv(first$metrics, f1m)
  write_tsv(first$pmetrics, f1pm)
  write_tsv(first$obs_as_query, f1oq)
  write_tsv(first$true_as_query, f1tq)
} else {
  first <- get_id_files(f1tq, f1oq, f1m, f2pm)
}

if (!file.exists(f2tq)) {
  second <- main("ptools", "2-Second_pass")
  write_tsv(second$metrics, f2m)
  write_tsv(second$pmetrics, f2pm)
  write_tsv(second$obs_as_query, f2oq)
  write_tsv(second$true_as_query, f2tq)
} else {
  second <- get_id_files(f2tq, f2oq, f2m, f2pm)
}

combined_tq <- bind_rows(
  mutate(first$true_as_query, pass = "First pass"),
  mutate(second$true_as_query, pass = "Second pass"),
)

combined_oq <- bind_rows(
  mutate(first$obs_as_query, pass = "First pass"),
  mutate(second$obs_as_query, pass = "Second pass"),
)

sources <- unique(combined_tq$source)

plot_helper <- function(tb, palette) {
  tb |> ggplot(aes(y = similarity, fill = source)) +
    geom_boxplot(size = 1) +
    scale_fill_paletteer_d(palette) +
    M$default_theme +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    facet_wrap(~pass)
}

GRAPHS$compared_tq <- plot_helper(combined_tq, "rcartocolor::Pastel")
attr(GRAPHS$compared_tq, "width") <- 18
GRAPHS$compared_oq <- plot_helper(combined_oq, "RColorBrewer::Set3")
attr(GRAPHS$compared_oq, "width") <- 18

compare_helper <- function(current, name, paired = TRUE) {
  query_sims <- lapply(sources, \(x) {
    filter(current, source == x)$similarity
  }) %>% `names<-`(sources)
  if (paired) {
    fn2 <- \(x, y) wilcox.test(x, y, alternative = "greater", paired = TRUE)
  } else {
    fn2 <- \(x, y) wilcox.test(x, y, alternative = "greater")
  }
  tests <- bind_rows(
    test_all_pairs(query_sims, \(x, y) wilcox.test(x, y), two_sided = TRUE),
    test_all_pairs(query_sims, fn2, "greater")
  ) |> get_adjusted_p()
  conclusion <- conclude_one_sided(tests)
  format <- pairwise_conclusion2gt(conclusion)
  TABLES[[glue("{name}_all_tests")]] <<- conclusion
  TABLES[[glue("{name}_significant")]] <<- format
}

to_test <- list(
  ptools_true_as_query = combined_tq,
  ptools_obs_as_query = combined_oq
)
temp <- c(TRUE, FALSE)

lapply(seq_along(to_test), \(x) {
  cur <- names(to_test)[x]
  compare_helper(to_test[[cur]], cur, temp[x])
})


combined_m <- bind_rows(
  mutate(first$metrics, pass = "First"),
  mutate(second$metrics, pass = "Second"),
) |> mutate(tpr = round(true_positive_ratio, 2))


combined_pm <- bind_rows(
  mutate(first$pmetrics, pass = "First"),
  mutate(second$pmetrics, pass = "Second"),
) |> mutate(tpr = round(true_positive_ratio, 2))

graph_helper <- function(data, palette) {
  data |> ggplot(aes(
    y = n_true_positives, x = source,
    fill = source, pattern = pass
  )) +
    geom_bar_pattern(
      position = "dodge", stat = "identity", pattern_density = 0.1,
      pattern_spacing = 0.03
    ) +
    ylab("Number of true positives") +
    M$default_theme +
    geom_text(aes(label = tpr, y = n_true_positives + 15),
      position = position_dodge(width = 1), vjust = 0,
      size = 4
    ) +
    scale_fill_paletteer_d(palette) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    scale_pattern_manual(values = c(First = "none", Second = "stripe"))
}


if ("ggpattern" %in% as_tibble(installed.packages())$Package) {
  library("ggpattern")
  combined_both <- bind_rows(
    mutate(combined_m, type = "Peptides"),
    mutate(combined_pm, type = "Proteins")
  )
  GRAPHS$ptools_true_positive <- graph_helper(combined_m, "ggthemes::Green_Orange_Teal")
  GRAPHS$ptools_true_positive_prot <- graph_helper(combined_pm, "ggthemes::Red_Blue_Brown")
  GRAPHS$ptools_true_positive_prot
  GRAPHS$ptools_both <- graph_helper(combined_both, "ggthemes::Green_Orange_Teal") + facet_wrap(~type)
  attr(GRAPHS$ptools_true_positive, "width") <- 17
  attr(GRAPHS$ptools_both, "width") <- 20
  attr(GRAPHS$ptools_true_positive_prot, "width") <- 17
}


TABLES$ptools_metrics <- combined_m |>
  select(-tpr) |>
  arrange(source) |>
  gt() |>
  cols_label(
    source = "Source",
    n_true_positives = "Number of true positives",
    true_positive_ratio = "True positive ratio",
  )

TABLES$ptools_pmetrics <- combined_pm |>
  select(-tpr) |>
  arrange(source) |>
  gt() |>
  cols_label(
    source = "Source",
    n_true_positives = "Number of true positives",
    true_positive_ratio = "True positive ratio",
  )

save(c(GRAPHS, TABLES), outdir)
