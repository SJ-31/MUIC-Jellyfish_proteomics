library(tidyverse)
library(glue)

# Sort proteins in the intersected searches output file into new groups


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
unify_groups <- function(combined_tib) {
  protein_record <- list()
  group_record <- list()
  counter <- 1
  for (index in seq_along(combined_tib$ProteinId)) {
    current <- combined_tib[index, ]
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
  return(protein_record)
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
  new <- lapply(colnames(currentID), collect_in_group, dataframe = currentID) %>%
    `names<-`(colnames(currentID))
  return(as_tibble(new))
}

resolve_records <- function(combined_tib, record) {
  all_groups <- names(record) %>%
    lapply(., group_frame,
      combined_frame = combined_tib,
      p_record = record
    ) %>%
    bind_rows() %>%
    mutate(num_proteins = unlist(lapply(ProteinId, num_prot))) %>%
    mutate(num_peps = unlist(lapply(peptideIds, num_prot))) %>%
    mutate(ProteinGroups = names(record))
  return(all_groups)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  combined_file <- args[1]
  output <- args[2]
  combined <- read.delim(combined_file, sep = "\t") %>%
    as_tibble()
  group_records <- unify_groups(combined)
  all_groups <- count_ids(combined, group_records)
  write_delim(all_groups, output, delim = "\t")
}
