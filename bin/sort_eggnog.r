library(seqinr)
library(tidyverse)
test <- TRUE
if (test) {
  unmatched <- "../results/jellyfish/1-First_pass/Unmatched/no_one_hits_degenerates_unmatched.tsv"
  ortho <- "../results/jellyfish/1-First_pass/Unmatched/eggNOG/no_one_hits_degenerates_unmatched.emapper.orthologs"
  hits <- "../results/jellyfish/1-First_pass/Unmatched/eggNOG/no_one_hits_degenerates_unmatched.emapper.hits"
  anno <- "../results/jellyfish/1-First_pass/Unmatched/eggNOG/no_one_hits_degenerates_unmatched.emapper.annotations"
  seed <- "../results/jellyfish/1-First_pass/Unmatched/eggNOG/no_one_hits_degenerates_unmatched.emapper.seed_orthologs"
  original_unmatched <- read_tsv(unmatched)
  ortho_df <- read_tsv(ortho, skip = 4)
  hit_df <- read_tsv(hits, skip = 4)
  anno_df <- read_tsv(anno, skip = 4)
  seed_df <- read_tsv(seed, skip = 5) %>% select(-c("evalue"))
  interpro <- "../results/jellyfish/1-First_pass/Unmatched/InterPro/no_one_hits_degenerates_unmatched-SCAN.tsv"
  interpro_df <- read_tsv(interpro)
}

clean_interpro <- function(df) {
  ip_df <- df %>%
    filter(evalue < 0.001) %>%
    ## filter(GO != "NA", GO != "-", pathways != "NA", pathways != ",") %>% # Too extreme?
    mutate(interpro_db = paste0(member_db, ": ", db_accession),
           interpro_db = recode(interpro_db,
                                "Coils: Coil" = "Coil",
                                "MobiDBLite: mobidb-lite" = "MobiDBLite"),
           interpro_description = paste0(description, ", ",
                                         interpro_description),
           interpro_evalue = evalue,
           GO = str_replace(GO, "NA", ""),
           pathways = str_replace(pathways, "NA", "")) %>%
    select(-c("sequence_md5", "length", "status", "date", "evalue",
              "db_accession", "description")) %>%
    group_by(query) %>%
    ungroup()
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
  if (num_searches == 1) return(list(OVERLAP1 = 1))
  while (m <= num_searches) {
    current <- query_df[m, ]
    while (n <= num_searches) {
      if (m == n && num_searches == 1) return(list(OVERLAP1 = 1))
      if (m == n) {
        n <- n + 1
        next
      }
      compare <- query_df[n, ]
      key <- paste0("OVERLAP", num_overlaps)
      if (compare$start <= current$stop) {
        overlapping[[key]] <- c(overlapping[[key]],m, n)
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
  current <- df[overlapping_indices, ]
  current <- distinct(current, interpro_db, .keep_all = TRUE) # interpro_db is
  # the column to consider for distinct
  predictive <- c("Coils", "MobiDBLite", "SUPERFAMILY", "PIRSF", "Gene3D",
                  "FunFam")
  curated <- c("Hamap", "PANTHER", "Pfam", "ProSiteProfiles", "CDD",
               "ProSitePatterns", "ProSiteProfiles", "SMART", "PRINTS",
               "SFLD", "NCBIfam")
  if (any(current$member_db %in% predictive) && any(current$member_db %in% curated)) {
    current <- current %>% filter(current$member_db %in% curated)
  }
  return(current)
}

distinct_from_csv <- function(csv_string) {
  split <- strsplit(csv_string, ",")
  split <- unique(unlist(split, use.names = FALSE))
  return(paste0(split[split != "-"], collapse = ","))
}

clean_annotations <- function(df) {
  # Clean up all overlapping annotations for the query sequences in "df"
  ip_df <- clean_interpro(df)
  queries <- ip_df$query %>% unique()
  cleaned <- lapply(queries, function(x) {
    current_query <- ip_df %>% filter(query == x) %>% arrange(start)
    current_overlaps <- get_overlaps(ip_df, x)
    current_sorted <- lapply(current_overlaps,
                             resolve_overlaps, df = current_query) %>%
      bind_rows() %>%
      distinct(., interpro_description, .keep_all = TRUE) %>%
      group_by(query) %>%
      mutate_at(., vars(-group_cols()), paste0, collapse = ",") %>%
      distinct() %>%
      mutate_at(., vars(-group_cols()), distinct_from_csv) %>%
      ungroup() %>%
      select(-c(start, stop))
    return(current_sorted)
  }) %>%
    bind_rows()
  return(cleaned)
}


interpro_cleaned <- clean_annotations(interpro_df)

# eggnog sorting
distinct_cols <- c("query", "evalue", "score", "bitscore",
                   "pident", "qcov", "scov")
joined <- inner_join(anno_df, seed_df,
                     by = join_by(x$`#query` == y$`#qseqid`)) %>%
  select(-c("qstart", "qend", "sstart", "send", "sseqid")) %>%
  as_tibble()
joined <- inner_join(interpro_cleaned, joined,
                     by = join_by(x$query == y$`#query`))

merged <- joined %>%
  lapply(., as.character) %>%
  as_tibble() %>%
  group_by(seed_ortholog) %>%
  mutate_at(., distinct_cols,
            paste0, collapse = ",") %>%
  distinct() %>%
  ungroup()
