library(tidyverse)
library(glue)

# Sort proteins in the intersected searches output file into new groups


check_group <- function(group_list, group_record) {
  group_list <- as.character(group_list)
  split <- strsplit(group_list, split = ",") %>% unlist()
  try <- intersect(names(group_record), split)
  if (length(try) > 0) {
    not_in <- setdiff(split, names(group_record))
    existing_group <- group_record[[names(group_record)[1]]]
    addition <- as.list(rep(existing_group, length(not_in))) %>%
      `names<-`(not_in)
    return(addition)
  }
  return(0)
}

# If a protein group has already been recorded, add the current protein and
#  the groups that haven't been recorded into the unified group
#   else, create a new unified group
unify_groups <- function(combined_tib) {
  protein_record <- list()
  counter <- 1
  for (index in seq_along(combined_tib$ProteinId)) {
    current <- combined_tib[index, ]
    current_groups <- strsplit((current$ProteinGroupId), ",") %>% unlist()
    current_groups <- current_groups[grep("NA", current_groups,
      invert = TRUE
    )]
    checked <- check_group(current_groups, protein_record)
    if (!is.list(checked)) {
      name <- glue("G{counter}")
      new_group <- rep(name, length(current_groups)) %>%
        `names<-`(current_groups)
      protein_record <- append(protein_record, new_group)
      counter <- counter + 1
    } else {
      protein_record <- append(protein_record, checked)
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

resolve_records <- function(combined_tib, records) {
  new <- combined_tib %>%
    mutate(Group = unlist(lapply(ProteinGroupId, function(x) {
      groups <- strsplit(x, ",")[[1]]
      found <- intersect(groups, names(records))[1]
      return(records[[found]])
    })))
  return(new)
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  combined_file <- args[1]
  output <- args[2]
  combined <- read.delim(combined_file, sep = "\t") %>%
    as_tibble()
  group_records <- unify_groups(combined)
  all_groups <- resolve_records(combined, group_records) %>%
    mutate(Identification_method = "standard_search")
  write_delim(all_groups, output, delim = "\t")
}
