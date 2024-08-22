library("glue")
library("ggpattern")
library("gridExtra")
library("grid")
library("paletteer")
library("tidyverse")
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  tools <- "/home/shannc/Bio_SDD/tools"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
  tools <- "/home/shannc/workflow/tools"
}
source(glue("{wd}/bin/R/analysis/all_analyses.r"))

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

correct_taxa_names <- function(x) {
  case_match(
    x,
    "Homo" ~ "Human",
    "Escherichia" ~ "E. coli",
    "Saccharomyces" ~ "Yeast",
  )
}

# Given sample distribution
# 65% w/w  Human
# 22.5% w/w Yeast
# 12.5% w/w E.coli
IDEAL_RATIO <- tibble(
  taxon = c("Human", "Yeast", "E. coli"),
  proportion = c(65, 22.5, 12.5),
  source = "spike-in ratio"
)

get_prop_compare <- function(tb, lfq, grouping_col = "GroupUP", tax = NULL) {
  if (!is.null(tax)) {
    tb <- tb |>
      inner_join(tax, by = join_by(ProteinId))
  }
  has_tax <- tb |>
    group_by(!!as.symbol(grouping_col)) |>
    summarise(
      Genus = nth(Genus, 1),
      check_genera = paste0(unique(Genus), collapse = ";")
    )
  genus_mismatch <- has_tax |> filter(grepl(";", check_genera))
  id_proportion <- local({
    t <- has_tax$Genus |> table()
    t / sum(t)
  })

  compare_props <- id_proportion |>
    table2tb(id_col = "taxon") |>
    mutate(
      proportion = n * 100, source = "identifications",
      taxon = correct_taxa_names(taxon)
    ) |>
    select(-n) |>
    bind_rows(IDEAL_RATIO)
  lfq <- merge_lfq(tb, "mean")
  if (!is.null(tax)) {
    lfq <- inner_join(lfq, select(tax, ProteinId, Genus), by = join_by(ProteinId))
  }
  lfq <- lfq |>
    inner_join(select(tb, ProteinId, !!as.symbol(grouping_col))) |>
    group_by(!!as.symbol(grouping_col)) |>
    summarise(
      log_intensity = mean(log_intensity, na.rm = TRUE),
      Genus = nth(Genus, 1)
    )

  lfq_genus <- lfq |>
    group_by(Genus) |>
    summarise(
      log_intensity = sum(log_intensity, na.rm = TRUE)
    ) |>
    mutate(taxon = correct_taxa_names(Genus))

  lfq_genus$log_intensity_ratio <- lfq_genus$log_intensity / sum(lfq_genus$log_intensity)

  compare_props <- compare_props %>% bind_rows(
    select(lfq_genus, log_intensity_ratio, taxon) |>
      rename(proportion = log_intensity_ratio) |>
      mutate(
        source = "intensity",
        proportion = proportion * 100
      )
  )
  compare_props
}


GRAPHS <- list()
outdir <- glue("{wd}/docs/figures/benchmark")
args <- list(
  prefix = "benchmark",
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  chosen_pass = "2-Second_pass"
)
reticulate::source_python(glue("{args$python_source}/helpers.py"))
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))

seq_map_path <- glue("{outdir}/bm_seq_map.tsv")

if (!file.exists(seq_map_path)) {
  SEQ_MAP <- read_tsv(glue("{wd}/results/benchmark/Databases/seq-header_mappings.tsv")) |>
    filter(!grepl("rev_", id)) |>
    mutate(organism = NA) |>
    get_organism("Unknown")
  write_tsv(SEQ_MAP, seq_map_path)
} else {
  SEQ_MAP <- read_tsv(seq_map_path)
}

cov_tests <- tibble()
main <- function(args, palette) {
  benchmark_dir <- glue("{wd}/results/{args$prefix}")
  run <- get_run(args$prefix, benchmark_dir)
  merged <- inner_join(run$first, run$second, by = join_by(header))
  test <- wilcox.test(merged$pcoverage_align.y, merged$pcoverage_align.x, alternative = "greater", paired = TRUE) |>
    htest2tb(data.name = glue("second greater, {args$prefix}"))
  cov_tests <<- bind_rows(cov_tests, test)

  lfq <- read_tsv(glue("{benchmark_dir}/{args$chosen_pass}/lfq_all.tsv"))
  data <- read_tsv(glue("{benchmark_dir}/{args$chosen_pass}/{args$prefix}_all_wcoverage.tsv")) |>
    filter(grepl("P", ProteinId)) |>
    group_by_unique_peptides() |>
    inner_join(lfq, by = join_by(ProteinId))

  tax <- read_tsv(glue("{benchmark_dir}/{args$chosen_pass}/{args$prefix}_taxonomy.tsv")) |>
    filter(ProteinId %in% data$ProteinId)

  # blast <- read_tsv(glue("{benchmark_dir}/{args$chosen_pass}/Unmatched/BLAST/accepted_queries.tsv")) |>
  #   select(ProteinId, MatchedPeptideIds)
  # blast <- blast %>% separate_longer_delim(MatchedPeptideIds, ";")
  # perc_file <- glue("{benchmark_dir}/{args$chosen_pass}/percolator_all.tsv")
  # perc <- read_tsv(perc_file) |>
  #   mutate(join_key = paste0(ProteinId, engine, collapse = "-"))

  # if (!"Genus" %in% colnames(perc)) {
  #   perc <- inner_join(perc, SEQ_MAP, by = join_by(x$ProteinId == y$id)) |>
  #     mutate(Genus = map_chr(organism, \(x) str_split_1(x, " ")[1]))
  # }

  # browser()
  # matching <- perc |>
  #   inner_join(blast, by = join_by(x$ProteinId == y$MatchedPeptideIds)) |>
  #   group_by(ProteinId.y, engine) |>
  #   summarise(
  #     peptideIds = paste0(peptideIds, collapse = ";"),
  #     `q-value` = median(`q-value`)
  #   ) |>
  #   mutate(join_key = paste0(ProteinId.y, engine, collapse = "-")) |>
  #   inner_join(SEQ_MAP, by = join_by(x$ProteinId.y == y$id))

  # # browser()
  # perc <- bind_rows(
  #   perc,
  #   filter(matching, !ProteinId.y %in% perc$ProteinId)
  # ) |> filter(!grepl("D", ProteinId))

  # left_join(perc, matching, by = join_by(join_key))


  ## Check quantification

  # Will need to consult about the best way to group_by and discard proteins that could be
  # in the same group
  # Multiple options
  # 1. Consider only unique groups,
  # 2. Consider only unique combinations of percolator groups CONFIRMED
  # 3. Consider proteins with unique peptide ids
  compare_props <- get_prop_compare(data, lfq, tax = tax)


  new_labels <- c(
    "Identification ratio" = "identifications",
    "Intensity ratio" = "intensity",
    "Spike-in ratio (true ratio)" = "spike-in ratio"
  )

  GRAPHS$prop_comparison <- compare_props |>
    mutate(proportion = round(proportion, 2)) |>
    ggplot(
      aes(source, proportion, color = source, fill = taxon)
    ) +
    geom_col(
      linewidth = 3, show.legend = c(fill = TRUE, color = FALSE),
    ) +
    scale_fill_paletteer_d(palette) +
    scale_color_manual(values = c(
      "identifications" = CATPUCCIN_LATTE$mauve,
      "intensity" = CATPUCCIN_LATTE$mauve,
      "spike-in ratio" = CATPUCCIN_LATTE$text
    )) +
    geom_text(aes(label = proportion),
      position = position_stack(vjust = 0.5),
      size = 5, color = "black"
    ) +
    ylab("Ratio (%)") +
    xlab("Source") +
    guides(fill = guide_legend(title = "Taxon")) +
    scale_x_discrete(
      labels = get_label_replacement(new_labels),
      limits = c("identifications", "spike-in ratio", "intensity")
    ) +
    theme(axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 15))

  ggsave(
    filename = glue("{wd}/docs/{args$prefix}_prop_comparison.png"),
    plot = GRAPHS$prop_comparison
  )
  list(graph = GRAPHS$prop_comparison, data = compare_props)
}

default <- main(args, "PNWColors::Shuksan2")
default_g <- default$graph + guides(fill = "none")

args$prefix <- "benchmark.calibrated"
calib <- main(args, "khroma::pale")

calib_g <- calib$graph + theme(axis.title.y = element_blank()) + guides(fill = "none")
args$prefix <- "benchmark.msconvert"

msconvert <- main(args, "khroma::pale")
msconvert_g <- msconvert$graph + theme(axis.title.y = element_blank())

all <- bind_rows(
  mutate(default$data, param = "default"),
  mutate(calib$data, param = "Calibrated"),
  mutate(msconvert$data, param = "msConvert"),
) |>
  mutate(proportion = round(proportion, 2)) |>
  filter(source != "spike-in ratio")

pos <- list(
  yeast = 12.5 / 2,
  human = (65 / 2) + 12.5,
  ecoli = (22.5 / 2) + 12.5 + 65
)
colors <- list(
  yeast = CATPUCCIN_LATTE$green,
  human = CATPUCCIN_LATTE$yellow,
  ecoli = CATPUCCIN_LATTE$red
)

midpoint <- mean(seq_along(unique(all$param)))
ref_size <- 30

ref_legend <- grid::legendGrob(
  c("Yeast (0.225)", "Human (0.65)", "E. coli (0.125)"),
  gp = grid::gpar(col = c(colors$yeast, colors$human, colors$ecoli)),
  vgap = 0.4,
  pch = 15
)

GRAPHS$all_lfq <- all |> ggplot(aes(param, proportion, fill = taxon)) +
  geom_col() +
  scale_fill_paletteer_d("PNWColors::Sailboat") +
  geom_text(aes(label = proportion),
    position = position_stack(vjust = 0.5),
    size = 5, color = "black"
  ) +
  annotate("segment",
    x = 0, y = 0, yend = 22.5,
    color = colors$yeast, size = ref_size,
  ) +
  annotate("segment",
    x = 0, y = 22.5, yend = 65 + 22.5,
    color = colors$human, size = ref_size,
  ) +
  annotate("segment",
    x = 0, y = 65 + 22.5, yend = 100,
    color = colors$ecoli, size = ref_size,
  ) +
  guides(custom = ggplot2::guide_custom(ref_legend, title = "True ratios")) +
  xlab("Run parameter") +
  ylab("Proportion") +
  theme(
    legend.title = element_text(face = "bold", size = 13),
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(size = 13),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
  ) +
  facet_wrap(~source)
attr(GRAPHS$all_lfq, "width") <- 12

# DO this for individual engines BUT maybe latter

join_denovo <- function(perc, blast) {
  matching <- perc |>
    inner_join(blast, by = join_by(x$proteinid == y$matchedpeptideids)) |>
    group_by(proteinid.y, engine) |>
    summarise(
      peptideids = paste0(peptideids, collapse = ";"),
      `q-value` = median(`q-value`)
    ) |>
    mutate(join_key = paste0(proteinid.y, engine, collapse = "-")) |>
    inner_join(seq_map, by = join_by(x$proteinid.y == y$id))

  perc <- bind_rows(
    perc,
    filter(matching, !proteinid.y %in% perc$proteinid)
  ) |> filter(!grepl("D", proteinid))
  perc <- left_join(perc, matching, by = join_by(join_key))
}

save(GRAPHS, outdir)
