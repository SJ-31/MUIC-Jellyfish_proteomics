library("tidyverse")
library("glue")


main <- function(args) {
  percolator_protein_files <- list.files(args$percolator_dir,
    pattern = "_percolator_proteins.tsv", full.names = TRUE
  )
  combined_results <- read_tsv(args$input)

  kept_seqs <- combined_results$unique_peptides |> flatten_by(";")
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

  protein_tbs <- local({
    get_engine <- function(file) {
      tb <- separate_longer_delim(read_tsv(file[[1]]), "ProteinId", ",") %>% mutate(
        engine = names(file),
        modifiedPeptideIds = peptideIds,
        peptideIds = map_chr(peptideIds, clean_peptide)
      )
      filtered <- dplyr::filter(tb, ProteinId %in% all_found_protein_ids) %>%
        mutate(num_peps = map_dbl(peptideIds, \(x) length(str_split_1(x, " "))))
      return(tb)
    }
    lmap(percolator_files, get_engine) %>% `names<-`(ENGINES)
  })

  seq_map <- get_seq_map(
    args$seq_map_path, all_found_protein_ids,
    kept_seqs
  )

  percolator_tb <- dplyr::bind_rows(protein_tbs)
  peptide_tb <- percolator_tb |>
    mutate(peptideIds = map_chr(peptideIds, clean_peptide)) |>
    separate_longer_delim(peptideIds, " ") |>
    distinct(engine, peptideIds, .keep_all = TRUE)
  dlfq <- read_tsv(args$dlfq_path) |>
    rename(ProteinId = protein, peptideIds = ion) |>
    select(ProteinId, peptideIds) |>
    filter(!peptideIds %in% peptide_tb$peptideIds) |>
    separate_longer_delim(ProteinId, ";")
  peptide_tb <- bind_rows(peptide_tb, dlfq)
  result <- list(
    percolator_all = percolator_tb,
    seq_map = seq_map,
    peptides = peptide_tb
  )
}


clean_peptide <- function(str) {
  str |>
    str_remove_all("\\[[0-9]+\\.[0-9]+\\]") |>
    str_remove_all("\\[[a-z \\:A-Z]+\\]")
}

get_seq_map <- function(seq_map_path, kept_ids, kept_seqs) {
  seq_map <- read_tsv(seq_map_path) |>
    filter(!grepl("rev_", id)) |>
    filter(id %in% kept_ids | seq %in% kept_seqs)
  seq_map
}


if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-p", "--percolator_dir"), type = "character", action = "store")
  parser <- add_option(parser, c("-o", "--outdir"), type = "character", action = "store")
  parser <- add_option(parser, c("-r", "--r_source"), type = "character", action = "store")
  parser <- add_option(parser, c("-i", "--input"), type = "character", action = "store")
  parser <- add_option(parser, c("-s", "--seq_map_path"), type = "character", action = "store")
  parser <- add_option(parser, c("-d", "--dlfq_path"), type = "character", action = "store")
  parser <- add_option(parser, c("-f", "--fdr"),
    type = "numeric",
    action = "store", default = 0.05
  )
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  results <- main(args)
  write_tsv(results$percolator_all, glue("{args$outdir}/percolator_all.tsv"))
  write_tsv(results$seq_map, glue("{args$outdir}/seq-header_map_found.tsv"))
  write_tsv(results$peptides, glue("{args$outdir}/percolator_peptide_map.tsv"))
}
