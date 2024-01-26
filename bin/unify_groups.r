library(tidyverse)
library(glue)

# Sort proteins in the intersected searches output file into new groups

check_group <- function(group_list, group_record) {
  group_list <- as.character(group_list)
  split <- strsplit(group_list, split = ";") %>% unlist()
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
    current_groups <- strsplit((current$ProteinGroupId), ";") %>% unlist()
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
  return(paste(dataframe[[group_name]], collapse = ";"))
}

num_prot <- function(prot_list) {
  return(
    prot_list %>% strsplit(";") %>% unlist() %>% length()
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
      groups <- strsplit(x, ";")[[1]]
      found <- intersect(groups, names(records))[1]
      return(records[[found]])
    })))
  return(new)
}

main <- function(search_type, group_prefix, input) {
  combined <- read.delim(input, sep = "\t") %>%
    as_tibble()
  group_records <- unify_groups(combined)
  all_groups <- resolve_records(combined, group_records) %>%
    mutate(ID_method = search_type)
  return(all_groups)
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-o", "--output"),
    type = "character",
    help = "Output file name"
  )
  parser <- add_option(parser, c("-i", "--input"),
    type = "character",
    help = "Input file name"
  )
  parser <- add_option(parser, c("-s", "--search_type"), type = "character")
  parser <- add_option(parser, c("-p", "--group_prefix"), type = "character")
  args <- parse_args(parser)
  m <- main(args$search_type, args$group_prefix, args$input)
  write_delim(m, args$output, delim = "\t")
}
