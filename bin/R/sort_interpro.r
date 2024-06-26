library(seqinr)
library(tidyverse)
library(glue)

clean_interpro <- function(df) {
  # Reformat names for interpro member dbs and remove those that are not valid
  ip_df <- df %>%
    filter(evalue < 0.001) %>%
    ## filter(GO != "NA", GO != "-", pathways != "NA", pathways != ",") %>% # Too extreme?
    mutate(
      interpro_db = paste0(member_db, ": ", db_accession),
      interpro_db = recode(interpro_db,
                           "Coils: Coil" = "Coil",
                           "MobiDBLite: mobidb-lite" = "MobiDBLite"
      ),
      interpro_description = paste0(
        description, "; ",
        interpro_description
      ),
      interpro_evalue = evalue,
      GO = str_replace(GO, "NA", ""),
      GO = str_replace_all(GO, "\\|", ";"),
      pathways = str_replace(pathways, "NA", "")
    ) %>%
    select(-c(
      "sequence_md5", "length", "status", "date", "evalue",
      "db_accession", "description"
    )) %>%
    group_by(query) %>%
    ungroup() %>%
    mutate(inferred_by = "interpro") %>%
    rename(., all_of(c("interpro_pathways" = "pathways"))) %>%
    select(-c("member_db"))
  return(ip_df)
}

get <- function(list, query, default) {
  return(list[[query]] %||% default)
}

# How to decide between predictions?
# Could prioritize certain DBs, high score, low evalue
#
#
get_overlaps <- function(df, query_id) {
  # Given a dataframe with rows as annotations, return a list sorting
  # which rows have overlapping annotations
  overlapping <- list()
  query_df <- df %>%
    filter(query == query_id) %>%
    arrange(start)
  num_searches <- dim(query_df)[1]
  num_overlaps <- 1
  m <- 1
  n <- 1
  if (num_searches == 1) {
    return(list(OVERLAP1 = 1))
  }
  while (m <= num_searches) {
    current <- query_df[m,]
    while (n <= num_searches) {
      if (m == n && num_searches == 1) {
        return(list(OVERLAP1 = 1))
      }
      if (m == n) {
        n <- n + 1
        next
      }
      compare <- query_df[n,]
      key <- paste0("OVERLAP", num_overlaps)
      if (compare$start <= current$stop) {
        overlapping[[key]] <- c(overlapping[[key]], m, n)
      } else if (compare$start > current$stop) {
        num_overlaps <- num_overlaps + 1
        key <- paste0("OVERLAP", num_overlaps)
        overlapping[[key]] <- c(overlapping[[key]], n)
        m <- n
        break
      }
      n <- n + 1
    }
    m <- m + 1
  }
  return(overlapping)
}

resolve_overlaps <- function(overlapping_indices, df) {
  # Resolve the overlapping annotations given by the specified rows
  # Intended for use with output of
  # "get_overlaps" function
  current <- df[overlapping_indices,]
  current <- distinct(current, interpro_db, .keep_all = TRUE) # interpro_db is
  # the column to consider for distinct
  predictive <- c(
    "Coils", "MobiDBLite", "SUPERFAMILY", "PIRSF", "Gene3D",
    "FunFam"
  )
  curated <- c(
    "Hamap", "PANTHER", "Pfam", "ProSiteProfiles", "CDD",
    "ProSitePatterns", "ProSiteProfiles", "SMART", "PRINTS",
    "SFLD", "NCBIfam"
  )
  if (any(current$member_db %in% predictive) && any(current$member_db %in% curated)) {
    current <- current %>% filter(current$member_db %in% curated)
  }
  return(current)
}

distinct_from_csv <- function(csv_string) {
  split <- strsplit(csv_string, ";")
  split <- unique(unlist(split, use.names = FALSE))
  return(paste0(split[split != "-"], collapse = ";"))
}

clean_annotations <- function(df) {
  # Clean up all overlapping annotations for the query sequences in "df"
  ip_df <- clean_interpro(df)
  queries <- ip_df$query %>% unique()
  cleaned <- lapply(queries, function(x) {
    current_query <- ip_df %>%
      filter(query == x) %>%
      arrange(start)
    current_overlaps <- get_overlaps(ip_df, x)
    current_sorted <- lapply(current_overlaps,
                             resolve_overlaps,
                             df = current_query
    ) %>%
      bind_rows() %>%
      distinct(., interpro_description, .keep_all = TRUE) %>%
      group_by(query) %>%
      mutate_at(., vars(-group_cols()), paste0, collapse = ";") %>%
      distinct() %>%
      mutate_at(., vars(-group_cols()), distinct_from_csv) %>%
      ungroup() %>%
      select(-c(start, stop, interpro_evalue))
    return(current_sorted)
  }) %>%
    bind_rows()
  return(cleaned)
}

main <- function(args) {
  interpro_df <- read_tsv(args$interpro_results)
  eggnog_df <- read_tsv(args$eggnog_unmatched)
  # Note that the true sequence of the protein that
  # the peptide has mapped to is unknown
  cleaned <- clean_annotations(interpro_df)
  joined <- inner_join(select(eggnog_df, -inferred_by), cleaned,
                       by = join_by(x$ProteinId == y$query)
  )
  blanks_removed <- unlist(lapply(joined$interpro_accession, function(x) {
    if (x == "") {
      return(NA)
    } else {
      return(x)
    }
  }))
  joined$interpro_accession <- blanks_removed
  still_unmatched <- eggnog_df[!(
    eggnog_df$ProteinId %in% joined$ProteinId),] %>%
    select(-inferred_by)
  if (purrr::pluck_exists(joined, "query")) {
    joined <- dplyr::select(joined, -query)
  }
  return(list(unmatched = still_unmatched, all = joined))
}


writeUnmatched <- function(unmatched_peptides, unmatched, output) {
  # Find unmatched peptides
  unmatched_peptides <- filter(unmatched_peptides, ProteinId %in% unmatched$ProteinId) %>%
    rename(seq = peptideIds)
  unmatched <- filter(unmatched, !ProteinId %in% unmatched_peptides$ProteinId) %>%
    select(c(ProteinId, seq, header))
  to_write <- bind_rows(unmatched_peptides, unmatched) %>%
    mutate(header = paste0(ProteinId, " ", header))
  write.fasta(
    sequences = as.list(to_write$seq),
    names = to_write$header,
    file.out = output,
    nbchar = 100,
    open = "w", as.string = TRUE
  )
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--interpro_results"), type = "character")
  parser <- add_option(parser, c("-o", "--output"), type = "character")
  parser <- add_option(parser, c("-p", "--unmatched_peptides"), type = "character")
  parser <- add_option(parser, c("-u", "--eggnog_unmatched"), type = "character")
  parser <- add_option(parser, c("-f", "--final_unmatched"), type = "character")
  parser <- add_option(parser, c("-a", "--final_unmatched_fasta"), type = "character")
  args <- parse_args(parser)
  results <- main(args)
  unmatched_peptides <- read_tsv(args$unmatched_peptides)
  write_tsv(results$unmatched, args$final_unmatched)
  writeUnmatched(unmatched_peptides, results$unmatched,
                 args$final_unmatched_fasta
  )
  write_tsv(results$all, args$output)
}
