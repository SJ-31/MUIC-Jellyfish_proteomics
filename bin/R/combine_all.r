library(seqinr)
library(reticulate)
library(tidyverse)
library(glue)

OTHER_EGGNOG_COLS <- c(
  "BRITE", "CAZy", "BiGG_Reaction", "seed_ortholog",
  "eggNOG_OGs", "COG_category", "Description", "Preferred_name",
  "EC", "PFAMs"
)

FIRST_COLS <- c(
  "ProteinId", "ProteinGroupId", "header",
  "GroupUP", "MatchedPeptideIds", "NCBI_ID", "UniProtKB_ID",
  "organism", "lineage", "ID_method", "inferred_by"
)

#' Count of number unique GO terms assigned to each protein, and
#' determine the maximum semantic value of the GO terms of that protein
#'
count_go <- function(go_vector) {
  return(lapply(go_vector, \(x) {
    row <- tibble(GO_count = 0)
    if (!is.na(x)) {
      split <- unique(str_split_1(x, ";"))
      row$GO_count <- length(split)
    }
    return(row)
  }) %>% bind_rows())
}

#' Correct the "NCBI_ID" entry for database proteins that were originally derived
#' from UniProt
correct_ids <- function(tb) {
  from_uniprot <- local({
    mapped <- map2_lgl(
      # Due to how the NCBI_ID column was generated, proteins from
      # UniProt would have the same entry in these two columns
      tb$UniProtKB_ID, tb$NCBI_ID,
      \(x, y) str_detect(x, y)
    )
    replace_na(mapped, FALSE)
  })
  tb[from_uniprot, ]$NCBI_ID <- NA
  return(tb)
}

##'   Remove description from a vector of GO terms separated by ";"
clean_go_str <- function(go_vector) {
  go_vector %>%
    lapply(., \(x) {
      if (is.na(x)) {
        return(x)
      }
      s <- str_split_1(x, pattern = ";|,") %>% keep(grepl("GO:", .))
      r <- unlist(lapply(s, gsub, pattern = "_.*", replacement = "")) %>% unique()
      return(paste0(r, collapse = ";"))
    }) %>%
    unlist()
}


KEEP_AS_CHAR <- c(
  "ProteinId", "header", "NCBI_ID", "UniProtKB_ID", "organism", "ProteinGroupId",
  "lineage", "GO", "GO_evidence", "KEGG_Genes", "PANTHER", "MatchedPeptideIds",
  "peptideIds", "ID_method", "inferred_by", "seq",
  "seed_ortholog", "eggNOG_OGs", "COG_category", "Description",
  "Preferred_name", "EC", "KEGG_ko", "PFAMs", "KEGG_Pathway",
  "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
  "CAZy", "BiGG_Reaction", "interpro_accession",
  "interpro_description", "interpro_pathways", "interpro_db"
)

ANNO_COLS <- KEEP_AS_CHAR %>% discard(., \(x) {
  x %in% c("seq", "peptideIds", "ID_method", "inferred_by")
})

sort_vals <- function(values) {
  if (is.na(values)) {
    return(1)
  }
  v <- str_split_1(values, ";|,") %>% as.double()
  if (length(v) > 1) {
    v <- sort(v)[2]
  }
  return(v)
}

load_file <- function(path) {
  df <- read_tsv(path)
  to_character <- colnames(df)[!(colnames(df) %in% KEEP_AS_CHAR)]
  df <- df %>%
    mutate_at(., to_character, as.character) %>%
    mutate(across(where(\(x) all(is.na(x))), as.character))
  return(df)
}

#' Parse GO evidence codes for all entries that have them
#' GOs added by eggNOG and interpro are automatically labelled IEA
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

checkMatchedPeps <- function(tb) {
  return(tb$MatchedPeptideIds %>%
    map_lgl(\(x) grepl("P", x)) %>%
    any())
}


LOGFILE <- "./combine_all.log"

main <- function(args) {
  downloads <- load_file(args$download) %>%
    mutate(header = replace(header, header == "NaN", "unknown"))
  if (!file.size(args$eggnog) == 0) {
    eggnog <- load_file(args$eggnog)
    interpro <- load_file(args$interpro)
    combined <- bind_rows(downloads, eggnog, interpro)
  } else {
    combined <- downloads
  }

  combined[combined == ""] <- NA
  if ("NCBI_ID" %in% colnames(combined) && "UniProtKB_ID" %in% colnames(combined)) {
    combined <- correct_ids(combined)
  }
  combined <- combined %>% mutate(
    num_peps = str_count(peptideIds, ";") + 1,
    mass = map_dbl(mass, function(x) {
      if (x == "-" | is.na(x)) {
        return(NA)
      } else {
        return(as.double(x))
      }
    })
  )

  combined <- combined %>%
    mutate_all(~ replace(., . == "-", NA)) %>%
    mutate_all(~ replace(., . == NaN, NA)) %>%
    mutate(
      length = as.double(gsub("unknown", NA, length)),
      unique_peptides = purrr::map_chr(peptideIds, \(x) {
        if (is.na(x)) {
          return(NA)
        }
        x <- str_split_1(x, ";") %>% map_chr(clean_peptide)
        return(paste0(unique(x), collapse = ";"))
      }),
      num_unique_peps = str_count(unique_peptides, ";") + 1,
      ID_method = purrr::map_chr(
        ID_method,
        \(x) ifelse(grepl(";", x), "both", x)
      )
    )

  ## Filter by FDR and PEP
  ## If a (standard) protein has at least two peptides are lower than the fdr threshold, it is kept
  combined <- combined %>%
    dplyr::mutate(
      q_adjust = map_dbl(`q.value`, sort_vals),
      pep_adjust = map_dbl(posterior_error_prob, sort_vals)
    )
  combined <- combined |>
    filter(q_adjust <= args$fdr) %>%
    filter(pep_adjust <= args$pep_thresh)

  redundant <- local({
    cols <- sapply(
      colnames(combined),
      function(x) all(is.na(combined[[x]]))
    )
    names(cols[cols])
  })

  combined <- select(combined, -all_of(redundant))

  ## Map pfam domains to GO and get GO IDs
  if ("GO" %in% colnames(combined)) {
    combined <- combined |> mutate(GO_evidence = apply(., 1, getEvidence))
    cat("BEGIN: mapping pfam domains to GO IDs\n", file = LOGFILE, append = TRUE)
    tryCatch(
      expr = {
        pfam_env <- new.env()
        source_python(glue("{args$python_source}/map2go.py"), envir = pfam_env)
        combined <- pfam_env$map_all_db(
          to_annotate = as.data.frame(combined),
          p2g_path = args$pfam2go,
          k2g_path = args$kegg2go,
          ec2g_path = args$ec2go,
          i2g_path = args$interpro2go,
          pfam_db_path = args$pfam_db
        ) %>%
          as_tibble() %>%
          dplyr::mutate(GO_IDs = clean_go_str(GO))
        gc <- count_go(combined$GO_IDs)
        combined$GO_counts <- gc$GO_count
      },
      error = \(cnd) {
        print(reticulate::py_last_error())
        stop("Caught reticulate error")
      }
    )
    cat("COMPLETE: mapping pfam domains to GO IDs\n", file = LOGFILE, append = TRUE)
  }
  if ("GO" %in% colnames(combined)) {
    combined$GO_slims <- combined$GO_IDs %>%
      lapply(
        .,
        \(x) {
          slims_from_go(x,
            go_path = args$go_path,
            go_slim_path = args$go_slim_path
          )
        }
      ) %>%
      as.character()
  }

  ## Map KEGG Genes to KEGG pathways
  if (any(grepl("KEGG", colnames(combined)))) {
    cat("BEGIN: mapping kegg genes to pathways\n", file = LOGFILE, append = TRUE)
    kegg_env <- new.env()
    source_python(glue("{args$python_source}/map2kegg.py"), envir = kegg_env)
    tryCatch(
      expr = {
        combined <- kegg_env$mapInDf(as.data.frame(combined), "pathway", "KEGG_Pathway")
        combined <- kegg_env$mapInDf(as.data.frame(combined), "module", "KEGG_Module")
        combined <- kegg_env$mapInDf(as.data.frame(combined), "enzyme", "EC")
        combined <- kegg_env$mapInDf(as.data.frame(combined), "ko", "KEGG_ko")
      },
      error = \(cnd) {
        print(reticulate::py_last_error())
        stop("Caught reticulate error")
      }
    )
    cat("COMPLETE: mapping kegg genes to pathways\n", file = LOGFILE, append = TRUE)
  }

  ## Arrange columns
  rename_lookup <- c(
    eggNOG_preferred_name = "Preferred_name",
    eggNOG_description = "Description"
  )
  if (!"organism" %in% colnames(combined)) {
    combined$organism <- NA
  }
  combined <- get_organism(combined, denovo_org = str_replace_all(args$denovo_org, "_", " "))
  combined <- combined %>%
    dplyr::rename(any_of(rename_lookup)) %>%
    relocate(where(is.numeric),
      .after = where(is.character)
    ) %>%
    relocate(c("q.value", "posterior_error_prob"), .before = "q_adjust") %>%
    relocate(any_of(c("peptideIds", "seq")), .after = where(is.numeric)) %>%
    relocate(any_of(contains("interpro"))) %>%
    relocate(any_of(contains("KEGG"))) %>%
    relocate(any_of(c(
      contains("eggNOG"),
      "seed_ortholog", "COG_category"
    ))) %>%
    relocate(any_of(c(contains("GO"), "PANTHER"))) %>%
    relocate(contains("is_blast"), .before = where(is.numeric)) %>%
    relocate(any_of(FIRST_COLS))
  combined <- group_by_unique_peptides(combined)
  return(combined)
}


if (sys.nframe() == 0) { # Won't run if the script is being sourced
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, "--output", type = "character")
  parser <- add_option(parser, "--interpro", type = "character")
  parser <- add_option(parser, "--eggnog", type = "character")
  parser <- add_option(parser, "--downloads", type = "character")
  parser <- add_option(parser, "--denovo_org", type = "character")
  parser <- add_option(parser, "--fdr", type = "double")
  parser <- add_option(parser, "--pep_thresh", type = "double")
  parser <- add_option(parser, "--interpro2go", type = "character")
  parser <- add_option(parser, "--kegg2go", type = "character")
  parser <- add_option(parser, "--ec2go", type = "character")
  parser <- add_option(parser, "--pfam2go", type = "character")
  parser <- add_option(parser, "--go_path", type = "character")
  parser <- add_option(parser, "--go_slim_path", type = "character")
  parser <- add_option(parser, "--pfam_db", type = "character")
  parser <- add_option(parser, "--r_source", type = "character")
  parser <- add_option(parser, "--python_source", type = "character")
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  source(glue("{args$r_source}/GO_helpers.r"))
  reticulate::source_python(glue("{args$python_source}/helpers.py"))
  results <- main(args)
  results <- distinct(results)
  write_tsv(results, args$output)
}
