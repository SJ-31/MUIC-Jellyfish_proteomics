library("glue")
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


data <- read_tsv(glue("{args$benchmark}/{CHOSEN_PASS}/benchmark_all_wcoverage.tsv"))
tax <- read_tsv(glue("{args$benchmark}/{CHOSEN_PASS}/benchmark_taxonomy.tsv"))
data_rep <- read_tsv(glue("{args$benchmark}/{CHOSEN_PASS}/benchmark_all_wrepresentatives.tsv"))

# Given sample distribution
# 65% w/w  Human
# 22.5% w/w Yeast
# 12.5% w/w E.coli

## Check quantification
has_tax <- inner_join(data_rep, tax, by = join_by(ProteinId))
id_proportion <- local({
  t <- has_tax$Genus |> table()
  t / sum(t)
})
id_proportion


# Will need to consult about the best way to group_by and discard proteins that could be
# in the same group
# Multiple options
# 1. Consider only unique groups,
# 2. Consider only unique combinations of percolator groups
# 3. Consider proteins with unique peptide ids
tb <- data
# tb <- tb |> distinct(ProteinGroupId, .keep_all = TRUE)
lfq <- mergeLfq(tb, "mean") |>
  inner_join(select(has_tax, ProteinId, Genus), by = join_by(ProteinId))
lfq <- lfq |>
  group_by(Genus) |>
  nest() |>
  mutate(
    log_intensity = map_dbl(data, \(x) sum(x$log_intensity, na.rm = TRUE)),
  )
lfq$log_intensity_ratio <- lfq$log_intensity / sum(lfq$log_intensity)
lfq
# Get representatives then do ratio again
