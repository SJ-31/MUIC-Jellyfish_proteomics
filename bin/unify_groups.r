library(tidyverse)
library(glue)

# Sort proteins in the intersected searches output file into new groups

args <- commandArgs(trailingOnly = TRUE)
combined_file <- args[1]
output <- args[2]

combined <- read.delim(combined_file, sep = "\t") %>%
  as_tibble()

new_groups <- list()
seq_along(combined)

check_group <- function(group_list, group_record) {
  group_list <- as.character(group_list)
  split <- strsplit(group_list, split = ",") %>% unlist()
  for (index in seq_along(group_record)) {
    try <- intersect(group_record[[index]], split)
    if (length(try) > 0) {
      return(names(group_record[index]))
    }
  }
  return(0)
}

# If a protein group has already been recorded, add the current protein and
#  the groups that haven't been recorded into the unified group
#   else, create a new unified group
protein_record <- list()
group_record <- list()
counter <- 1
for (index in seq_along(combined$ProteinId)) {
  current <- combined[index, ]
  current_groups <- strsplit((current$ProteinGroupId), ",") %>% unlist()
  current_groups <- current_groups[grep("NA", current_groups,
    invert = TRUE
  )]
  checked <- check_group(current_groups, group_record)
  if (checked == 0) {
    name <- glue("G{counter}")
    group_record[[name]] <- current_groups
    protein_record[[name]] <- current$ProteinId
    counter <- counter + 1
  } else {
    group_record[[checked]] <- append(
      group_record[[checked]],
      current_groups
    ) %>% unique()
    protein_record[[checked]] <- append(
      protein_record[[checked]],
      current$ProteinId
    )
  }
}

collect_in_group <- function(group_name, dataframe) {
  return(paste(dataframe[[group_name]], collapse = ","))
}

num_prot <- function(prot_list) {
  return(
    prot_list %>% strsplit(",") %>% unlist() %>% length()
  )
}

group_frame <- function(pId, combined_frame, p_record) {
  currentID <- combined_frame %>% filter(ProteinId %in% p_record[[pId]])
  combined_row <- data.frame(
    ProteinId = collect_in_group("ProteinId", currentID),
    q.value = collect_in_group("q.value", currentID),
    original_groups = collect_in_group("ProteinGroupId", currentID),
    posterior_error_prob = collect_in_group("posterior_error_prob", currentID),
    peptideIds = collect_in_group("peptideIds", currentID),
    header = collect_in_group("header", currentID)
  )
  return(combined_row %>% as_tibble())
}

all_groups <- names(protein_record) %>%
  lapply(., group_frame,
    combined_frame = combined,
    p_record = protein_record
  ) %>%
  bind_rows() %>%
  mutate(num_proteins = unlist(lapply(ProteinId, num_prot))) %>%
  mutate(num_peps = unlist(lapply(peptideIds, num_prot))) %>%
  mutate(ProteinGroups = names(protein_record))

write_delim(all_groups, output, delim = "\t")
