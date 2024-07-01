PALETTE <- "PNWColors::Shuksan2"
library("glue")
library("gridExtra")
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

args <- list(
  benchmark = glue("{wd}/results/benchmark"),
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo"),
  go_tm_dir = glue("{wd}/data/reference/.go_texts")
)
CHOSEN_PASS <- "1-First_pass"
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))


data <- read_tsv(glue("{args$benchmark}/{CHOSEN_PASS}/benchmark_all_wcoverage.tsv")) |> group_by_unique_peptides()
tax <- read_tsv(glue("{args$benchmark}/{CHOSEN_PASS}/benchmark_taxonomy.tsv"))

ideal_ratio <- tibble(
  taxon = c("Human", "Yeast", "E. coli"),
  proportion = c(65, 22.5, 12.5),
  source = "spike-in ratio"
)
# Given sample distribution
# 65% w/w  Human
# 22.5% w/w Yeast
# 12.5% w/w E.coli

## Check quantification
# has_tax <- inner_join(data_rep, tax, by = join_by(ProteinId))
has_tax <- data |>
  inner_join(tax, by = join_by(ProteinId)) |>
  group_by(GroupUP) |>
  summarise(
    Genus = nth(Genus, 1),
    check_genera = paste0(unique(Genus), collapse = ";")
  )
genus_mismatch <- has_tax |> filter(grepl(";", check_genera))

id_proportion <- local({
  t <- has_tax$Genus |> table()
  t / sum(t)
}) # Pretty close, NICE
# Another win for grouping by unique shared peptides

correct_taxa_names <- function(x) {
  case_match(
    x,
    "Homo" ~ "Human",
    "Escherichia" ~ "E. coli",
    "Saccharomyces" ~ "Yeast"
  )
}

compare_props <- id_proportion |>
  table2tb(id_col = "taxon") |>
  mutate(
    proportion = n * 100, source = "identifications",
    taxon = correct_taxa_names(taxon)
  ) |>
  select(-n) |>
  bind_rows(ideal_ratio)

# Will need to consult about the best way to group_by and discard proteins that could be
# in the same group
# Multiple options
# 1. Consider only unique groups,
# 2. Consider only unique combinations of percolator groups CONFIRMED
# 3. Consider proteins with unique peptide ids
# tb <- tb |> distinct(ProteinGroupId, .keep_all = TRUE)
lfq <- merge_lfq(data, "mean") |>
  inner_join(select(tax, ProteinId, Genus), by = join_by(ProteinId)) |>
  inner_join(select(data, ProteinId, GroupUP)) |>
  group_by(GroupUP) |>
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
  scale_fill_paletteer_d(PALETTE) +
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
  )


ggsave(
  filename = glue("{wd}/docs/benchmark_prop_comparison.png"),
  plot = GRAPHS$prop_comparison
)
