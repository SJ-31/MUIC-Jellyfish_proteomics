library("tidyverse")
library("glue")

main <- function(args) {
  percolator_protein_files <- list.files(args$percolator_dir,
    pattern = "_percolator_proteins.tsv", full.names = TRUE
  )
  combined_results <- read_tsv(args$input)

  denovo_transcriptome_ids <- combined_results$MatchedPeptideIds %>%
    discard(is.na) %>%
    lapply(., str_split_1, pattern = ";") %>%
    unlist()
  database_ids <- combined_results$ProteinId
  all_found_protein_ids <- c(database_ids, denovo_transcriptome_ids)

  ENGINES <- percolator_protein_files %>%
    map_chr(., \(x) {
      gsub("_percolator_proteins.tsv", "", x) %>%
        gsub(".*/", "", .)
    })

  names(percolator_protein_files) <- ENGINES
  percolator_files <- as.list(percolator_protein_files)

  engine_stats <- tibble()
  protein_tbs <- local({
    get_engine <- function(file) {
      tb <- tb_duplicate_at(read_tsv(file[[1]]), "ProteinId", ",") %>% mutate(engine = names(file))
      filtered <- dplyr::filter(tb, ProteinId %in% all_found_protein_ids) %>%
        mutate(num_peps = map_dbl(peptideIds, \(x) length(str_split_1(x, " "))))
      row <- tibble(
        engine = names(file),
        n_groups_all = length(unique(tb$ProteinGroupId)),
        n_groups_accepted = length(unique(filtered$ProteinGroupId)),
        n_proteins_all = nrow(tb),
        n_proteins_below_fdr = tb %>% filter(`q-value` <= args$fdr) %>% nrow(),
        n_proteins_accepted = nrow(filtered),
        # To be accepted, a protein needs to be matched by at least two other engines
        # And the q value of at least two psm matches needs to be below the fdr
        n_other_proteins = map_lgl(tb$ProteinId, \(x) str_detect(x, "D|T")) %>% sum(),
      ) %>%
        mutate(n_database = n_proteins_all - n_other_proteins)
      engine_stats <<- dplyr::bind_rows(engine_stats, row)
      return(filtered)
    }
    lmap(percolator_files, get_engine) %>% `names<-`(ENGINES)
  })

  seq_map <- get_seq_map(args$seq_map_path, args$unmatched_path, all_found_protein_ids)

  percolator_tb <- dplyr::bind_rows(protein_tbs)
  peptide_tb <- tb_duplicate_at(percolator_tb, "peptideIds",
    separator = " ", split_fn = split_peptides
  ) |> select(ProteinId, ProteinGroupId, peptideIds, engine)
  result <- list(
    percolator_all = percolator_tb, metrics = engine_stats,
    seq_map = seq_map,
    peptides = peptide_tb
  )
}


split_peptides <- function(peptide_str) {
  clean_peptide <- function(str) {
    str |>
      str_remove_all("\\[[0-9]+\\.[0-9]+\\]") |>
      str_remove_all("\\[[a-z \\:A-Z]+\\]")
  }
  clean_peptide(peptide_str) |> str_split_1(" ")
}

get_seq_map <- function(seq_map_path, unmatched_path, kept_ids) {
  seq_map <- read_tsv(seq_map_path) |> filter(id %in% kept_ids)
  unmatched <- read_tsv(unmatched_path) |>
    filter(ProteinId %in% kept_ids) |>
    rename(id = ProteinId, seq = peptideIds) |>
    mutate(seq = map_chr(seq, clean_peptide)) |>
    select(id, seq)
  bind_rows(seq_map, unmatched)
}


if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-p", "--percolator_dir"), type = "character", action = "store")
  parser <- add_option(parser, c("-o", "--outdir"), type = "character", action = "store")
  parser <- add_option(parser, c("-r", "--r_source"), type = "character", action = "store")
  parser <- add_option(parser, c("-i", "--input"), type = "character", action = "store")
  parser <- add_option(parser, c("-s", "--seq_map_path"), type = "character", action = "store")
  parser <- add_option(parser, c("-f", "--fdr"),
    type = "numeric",
    action = "store", default = 0.05
  )
  parser <- add_option(parser, c("-u", "--unmatched_path"), type = "character", action = "store")
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  results <- main(args)
  write_tsv(results$percolator_all, glue("{args$outdir}/percolator_all.tsv"))
  write_tsv(results$metrics, glue("{args$outdir}/engine_stats.tsv"))
  write_tsv(results$seq_map, glue("{args$outdir}/seq-header_map_found.tsv"))
  write_tsv(results$peptides, glue("{args$outdir}/percolator_peptide_map.tsv"))
}
