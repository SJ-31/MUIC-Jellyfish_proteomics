library(reticulate)
library(tidyverse)
library(glue)

OTHER_EGGNOG_COLS <- c(
  "BRITE", "CAZy", "BiGG_Reaction", "seed_ortholog",
  "eggNOG_OGs", "COG_category", "Description", "Preferred_name",
  "EC", "PFAMs"
)

FIRST_COLS <- c(
  "ProteinId", "ProteinGroupId", "header", "Group", "NCBI_ID",
  "UniProtKB_ID", "organism", "lineage", "ID_method", "Anno_method"
)
##
##
##


writeAlignments <- function(row, file_name) {
  header <- ifelse(row[["header"]] == "unknown", row[["ProteinId"]],
    row[["header"]]
  )
  write.fasta(row[["seq"]], header, open = "a", file.out = file_name)
  write.fasta(row[["alignment"]], "ALIGNED PEPTIDES",
    open = "a",
    file.out = file_name
  )
  cat("\n", file = file_name, append = TRUE)
}

cleanPeptide <- function(pep) {
  if (grepl("\\]|[a-z0-9.]|-", pep)) {
    pep <- str_to_upper(pep) %>%
      str_extract_all("[A-Z]") %>%
      unlist() %>%
      paste0(collapse = "")
  }
  return(pep)
}

organismFromHeader <- function(row) {
  organism <- row["organism"]
  if (is.na(organism)) {
    header <- row["header"]
    if (grepl("OS=", header)) {
      return(str_extract(header, "OS=([a-zA-Z]* [a-zA-Z]*) ", group = 1))
    }
    return(str_extract(header, "\\[([A-Z].*)\\]", group = 1))
  }
  return(organism)
}

KEEP_AS_CHAR <- c(
  "ProteinId", "header", "NCBI_ID", "UniProtKB_ID", "organism", "ProteinGroupId",
  "lineage", "GO", "GO_evidence", "KEGG_Genes", "PANTHER",
  "peptideIds", "ID_method", "Anno_method", "seq",
  "seed_ortholog", "eggNOG_OGs", "COG_category", "Description",
  "Preferred_name", "EC", "KEGG_ko", "PFAMs", "KEGG_Pathway",
  "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
  "CAZy", "BiGG_Reaction", "interpro_accession",
  "interpro_description", "interpro_pathways", "interpro_db"
)

ANNO_COLS <- KEEP_AS_CHAR %>% discard(., \(x) {
  x %in% c("seq", "peptideIds", "ID_method", "Anno_method")
})

sortVals <- function(values) {
  v <- str_split_1(values, ",") %>% as.double()
  if (length(v) > 1) {
    v <- sort(v)[2]
  }
  return(v)
}

loadFile <- function(path) {
  df <- read_tsv(path)
  to_character <- colnames(df)[!(colnames(df) %in% KEEP_AS_CHAR)]
  df <- df %>%
    mutate_at(., to_character, as.character) %>%
    mutate(across(where(\(x) all(is.na(x))), as.character))
  return(df)
}

getEvidence <- function(row) {
  if (is.na(row[["GO"]])) {
    return(NA)
  }
  if (!is.na(row[["GO_evidence"]])) {
    evidence <- row[["GO_evidence"]] %>%
      str_split_1(";") %>%
      sapply(gsub, pattern = ":.*", replacement = "", USE.NAMES = FALSE) %>%
      paste0(., collapse = ";")
    return(evidence)
  }
  return("IEA")
}

## goFromPanther <- function(df) {
##   # Use the PANTHER api to obtain GO terms for proteins
##   # that lack the latter but have the former
##   to_map <- filter(df, !is.na(PANTHER) & is.na(GO))
##   the_rest <- filter(df, !(ProteinId %in% to_map$ProteinId))
## }


main <- function(args) {
  downloads <- loadFile(args$download) %>%
    mutate(header = replace(header, header == "NaN", "unknown"))
  if (str_to_lower(args$is_denovo) == "true") {
    eggnog <- loadFile(args$eggnog)
    interpro <- loadFile(args$interpro)
    combined <- bind_rows(downloads, eggnog, interpro)
  } else {
    combined <- downloads
  }
  combined[combined == ""] <- NA
  combined <- combined %>% mutate(
    num_peps = sapply(peptideIds, function(x) {
      return(str_count(x, ",") + 1)
    }, USE.NAMES = FALSE),
    mass = sapply(mass, function(x) {
      if (x == "-" || is.na(x)) {
        return(NA)
      } else {
        return(as.double(x))
      }
    }, USE.NAMES = FALSE)
  )
  redundant <- sapply(colnames(combined), function(x) {
    all(is.na(combined[[x]]))
  })

  combined <- combined %>%
    select(-c(names(redundant[redundant]))) %>%
    mutate_all(~ replace(., . == "-", NA)) %>%
    mutate_all(~ replace(., . == NaN, NA)) %>%
    mutate(
      GO_evidence = apply(., 1, getEvidence),
      length = as.double(gsub("unknown", NA, length)),
      unique_peptides = lapply(peptideIds, \(x) {
        x <- str_split_1(x, ",") %>%
          lapply(., cleanPeptide) %>%
          unlist()
        return(paste0(unique(x), collapse = ","))
      }) %>% unlist(),
      num_unique_peps = sapply(unique_peptides, \(x) {
        return(str_count(x, ",") + 1)
      }, USE.NAMES = FALSE),
    )

  ## Filter by FDR and PEP
  ## If a (standard) protein has at least two peptides are lower than the fdr
  ## threshold, it is kept
  combined <- combined %>%
    mutate(
      q_adjust = unlist(lapply(`q.value`, sortVals)),
      pep_adjust = unlist(lapply(posterior_error_prob, sortVals))
    ) %>%
    filter(q_adjust <= args$fdr) %>%
    filter(pep_adjust <= args$pep_thresh)

  ## Calculate coverage
  if (args$coverage) {
    source(glue("{args$r_source}/protein_coverage.r"))
    combined <- coverageCalc(combined)
    apply(filter(combined, !is.na(seq)), 1, writeAlignments,
      file_name = args$alignment_file
    )
  }
  ## Record modifications
  if (args$sort_mods) {
    source(glue("{args$r_source}/sort_mods.r"))
    combined <- sortModsMain(combined, FALSE)
  }

  ## Calculate empai using python script
  if (args$empai) {
    source_python(glue("{args$r_source}/emPAI.py"))
    combined <- py$calculate_emPAI(
      df = as.data.frame(combined),
      m_range = list(360L, 1600L)
    ) %>% as_tibble()
  }

  ## Merge with quantification data
  combined <- left_join(combined, read_tsv(args$directlfq),
    by = join_by(x$ProteinId == y$ProteinId)
  )
  combined <- left_join(combined, read_tsv(args$flashlfq),
    by = join_by(x$ProteinId == y$ProteinId)
  )

  ## Arrange columns
  combined <- combined %>%
    mutate(organism = apply(., 1, organismFromHeader)) %>%
    dplyr::rename(
      eggNOG_preferred_name = Preferred_name,
      eggNOG_description = Description
    ) %>%
    relocate(where(is.numeric),
      .after = where(is.character)
    ) %>%
    relocate(c("q.value", "posterior_error_prob"), .before = "q_adjust") %>%
    relocate(c("peptideIds", "SO_seq", "seq"), .after = where(is.numeric)) %>%
    relocate(contains("interpro")) %>%
    relocate(contains("KEGG")) %>%
    relocate(c(contains("eggNOG"), "seed_ortholog", "COG_category")) %>%
    relocate(c(contains("GO"), "PANTHER")) %>%
    relocate(contains("is_blast"), .before = where(is.numeric)) %>%
    relocate(all_of(FIRST_COLS))

  found <- colSums(!is.na(combined)) %>%
    lapply(., \(x) {
      x / nrow(combined) * 100
    }) %>%
    as.matrix()

  return(list(all = combined, f = found))
}



if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("--output"), type = "character")
  parser <- add_option(parser, c("--interpro"), type = "character")
  parser <- add_option(parser, c("--eggnog"), type = "character")
  parser <- add_option(parser, c("--downloads"), type = "character")
  parser <- add_option(parser, c("--fdr"), type = "double")
  parser <- add_option(parser, c("--alignment_file"), type = "character")
  parser <- add_option(parser, c("--pep_thresh"), type = "double")
  parser <- add_option(parser, c("--is_denovo"), type = "character")
  parser <- add_option(parser, c("--directlfq"), type = "character")
  parser <- add_option(parser, c("--flashlfq"), type = "character")
  parser <- add_option(parser, c("--coverage"),
    type = "character",
    default = TRUE,
    action = "store_true"
  )
  parser <- add_option(parser, c("--empai"),
    type = "character",
    default = TRUE,
    action = "store_true"
  )
  parser <- add_option(parser, c("--sort_mods"),
    type = "character",
    default = TRUE,
    action = "store_true"
  )
  parser <- add_option(parser, c("--r_source"), type = "character")
  args <- parse_args(parser)
  results <- main(args)
  write_tsv(results$all, args$output)
}
