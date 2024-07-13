library("tidyverse")
library("glue")
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
M$r_source <- glue("{M$wd}/bin/R")
source(glue("{M$r_source}/helpers.r"))

# Create mapping files
search_dir <- glue("{M$wd}/data/reference/search_queries")
previous_toxins_file <- glue("{M$wd}/data/reference/toxin_curated_pfamID.csv")

files <- list.files(search_dir, full.names = TRUE) |> discard(\(x) str_detect(x, "SearchResults"))

map_prev_fn <- swap_w_map(c(
  "neurotoxin;3FTxs" = "neurotoxin",
  "snake venom metalloproteinase" = "metalloprotease",
  "c-type lectin" = "c_type_lectin",
  "venom nerve growth factor" = "nerve_growth_factor",
  "disintegrin" = "metalloprotease",
  "amino acid oxidase" = "amino_acid_oxidase"
))
previous <- read_csv(previous_toxins_file) |>
  rename(Accession = pfam_accession, Group = keywords, Name = name) |>
  mutate(
    Group = map_prev_fn(Group),
    Source = "PFAM"
  )

toxin_rename_fn <- swap_w_map(c(
  "porins" = "pore_forming",
  "metallopeptidase" = "metalloprotease",
  "metalloproteases" = "metalloprotease",
  "neurotoxins" = "neurotoxin",
  "crisp" = "CRISP"
))


# Porins are subset of pore-forming proteins
toxins <- files |>
  lapply(\(x) {
    group <- basename(x) |>
      str_remove_all("interpro_|panther_|pfam_") |>
      str_remove("\\.tsv")
    read_tsv(x) |> mutate(Group = toxin_rename_fn(group))
  }) |>
  bind_rows() |>
  rename(Source = `Source Database`)

all_toxins <- bind_rows(toxins, previous) |>
  select(Accession, Name, Source, Group) |>
  distinct(.keep_all = TRUE)



write_tsv(all_toxins, glue("{M$wd}/data/reference/toxin_groups.tsv"))
