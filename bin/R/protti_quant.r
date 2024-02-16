library(tidyverse)
library(glue)
library(protti)

pth <- "./results/jellyfish/1-First_pass/Quantify/Mapped_scans/"

mostAbund <- function(grouped_df) {
  p <- pick(base_peptide)$base_peptide
  up <- unique(p)
  tab <- tabulate(match(p, up))
  abund <- up[tab == max(tab)]
  current <- grouped_df[cur_group_rows(), ]
  choices <- filter(current, base_peptide %in% abund)
  if (nrow(choices) > 1) {
    return(sample_n(choices, 1))
  }
  return(choices)
}

mapped_scans <- list.files(pth) %>%
  lapply(., function(x) {
    read_tsv(glue("{pth}/{x}"))
  }) %>%
  bind_rows() %>%
  filter(!is.na(protein))

## TESTING
mapped_scans <- mapped_scans[1:10000, ]

splitDuplicates <- function(row) {
  splits <- str_split_1(row[["protein"]], ";")
  stacked <- as.list(row) %>%
    lapply(., function(x) {
      return(rep(x, length(splits)))
    }) %>%
    as_tibble() %>%
    mutate(protein = splits)
  return(stacked)
}

sortDuplicates <- function(table) {
  keep_char <- c(
    "file", "scan", "base_peptide", "peptide", "protein",
    "engine", "precursor_id"
  )
  duplicates <- table %>% filter(grepl(";", protein))
  table <- table %>% filter(!grepl(";", protein))
  duplicates <- apply(duplicates, 1, splitDuplicates) %>%
    bind_rows() %>%
    mutate(across(!keep_char, as.double))
  bound <- bind_rows(list(duplicates, table))
  return(bound)
}


preprocess <- mapped_scans %>%
  select(c(protein, file, base_peptide, precursorIntensity, scan)) %>%
  group_by(scan) %>%
  mutate(
    precursorIntensity = log2(precursorIntensity),
    precursor_id = sapply(seq_along(nrow(.)), function(x) {
      glue("Pre{x}")
    })
  ) %>%
  summarize(mode = mostAbund(.))
preprocess <- preprocess$mode %>%
  ungroup() %>%
  sortDuplicates()

quantified <- calculate_protein_abundance(preprocess,
  sample = file, protein_id = protein, precursor = precursor_id,
  peptide = base_peptide, intensity_log2 = precursorIntensity, method = "iq",
  for_plot = FALSE
)
