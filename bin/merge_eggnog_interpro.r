library(seqinr)
library(tidyverse)
library(glue)
TEST <- TRUE
if (TEST) {
  dir <- "../results/jellyfish/1-First_pass/Unmatched"
  unmatched <- glue("{dir}/eggNOG/no_one_hits_degenerates_unmatched.tsv")
  orth <- glue("{dir}/eggNOG/no_one_hits_degenerates_unmatched.emapper.orthologs")
  hits <- glue("{dir}/eggNOG/no_one_hits_degenerates_unmatched.emapper.hits")
  anno <- glue("{dir}/eggNOG/no_one_hits_degenerates_unmatched.emapper.annotations")
  seed <- "../results/jellyfish/1-First_pass/Unmatched/eggNOG/no_one_hits_degenerates_unmatched.emapper.seed_orthologs"
  original_unmatched <- read_tsv(unmatched)
  ortho_df <- read_tsv(ortho, skip = 4)
  hit_df <- read_tsv(hits, skip = 4)
  anno_df <- read_tsv(anno, skip = 4)
  seed_df <- read_tsv(seed, skip = 5) %>% select(-c("evalue"))
  interpro <- "../results/jellyfish/1-First_pass/Unmatched/InterPro/no_one_hits_degenerates_unmatched-SCAN.tsv"
  interpro_df <- read_tsv(interpro)
}

## interpro_cleaned <- clean_annotations(interpro_df)

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
