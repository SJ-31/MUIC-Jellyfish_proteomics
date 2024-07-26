MAP_PFAMS <- TRUE
MAP_KEGG <- TRUE
UNIFY <- FALSE
MERGE_QUANT <- TRUE
TEST <- FALSE
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

HEADER_QUERIES <- list(
  unknown = "unknown|uncharacterized|predicted|unnamed|hypothetical",
  venom_component = "toxin|porin",
  transport = "pump|transport|channel",
  membrane = "membrane",
  translation = "ribosome|elongation|initiation|translation",
  catalytic_activity = "ase",
  cytoskeleton = "actin|tubulin|cytoskele|dynein|kinesin|catenin",
  other = ""
)


#' Count of number unique GO terms assigned to each protein, and
#' determine the maximum semantic value of the GO terms of that protein
#'
count_go <- function(go_vector) {
  return(lapply(go_vector, \(x) {
    row <- tibble(GO_count = 0, GO_max_sv = 0)
    if (!is.na(x)) {
      split <- unique(str_split_1(x, ";"))
      row$GO_count <- length(split)
      svs <- purrr::map_dbl(split, getSV)
      row$GO_max_sv <- max(svs)
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


#' When there are multiple quantified peptides for a protein, take
#' the top 3 peptides and average their values, reporting that as the
#' abundance of the protein
#' @param quant_name the name of the program used to generate the
#' quantification data e.g. "directlfq"
meanTop3 <- function(tb, quant_name) {
  tb <- tb %>%
    group_by(ProteinId) %>%
    nest() %>%
    apply(., 1, \(x) {
      top_three <- slice_max(x[["data"]], across(contains(quant_name)), n = 3)
      means <- top_three %>%
        summarise(across(contains(quant_name), \(x) mean(x, na.rm = TRUE)))
      joined_ids <- paste0(x[["data"]]$MatchedPeptideIds, collapse = ";")
      id <- x[["ProteinId"]]
      return(top_three[1, ] %>%
        select(-contains(quant_name)) %>%
        mutate(.,
          MatchedPeptideIds = joined_ids,
          ProteinId = id
        ) %>%
        bind_cols(means))
    }) %>%
    bind_rows()
  return(tb)
}


# Merge annotation tibble with quantification tibble
# Entries that have been matched by multiple peptides
# peptides (de novo, transcriptome etc.) in blast and interpro are
# extracted and handled differently: the values are averaged
# Compute average and median values across samples
merge_with_quant <- function(main_tb, quant_tb, quant_name) {
  has_multiple <- main_tb %>% filter(!is.na(MatchedPeptideIds))
  if (nrow(has_multiple) != 0) {
    # Attempt to find quantification for every peptide that was matched to a given protein
    split_up <- separate_longer_delim(has_multiple, "MatchedPeptideIds", ";")
    joined <- left_join(split_up, quant_tb, by = join_by(x$MatchedPeptideIds == y$ProteinId))
    has_multiple <- meanTop3(joined, quant_name = quant_name)
  }
  full_proteins <- main_tb %>% filter(is.na(MatchedPeptideIds))
  full_proteins <- left_join(full_proteins, quant_tb,
    by = join_by(x$ProteinId == y$ProteinId)
  )
  bound <- bind_rows(full_proteins, has_multiple)
  calcs <- bound %>%
    dplyr::rowwise() %>%
    dplyr::reframe(
      "{quant_name}_mean" := mean(c_across(contains(quant_name)), na.rm = TRUE),
      "{quant_name}_median" := median(c_across(contains(quant_name)), na.rm = TRUE)
    )
  return(bind_cols(bound, calcs))
}

unify_groups <- function(tb) {
  source_python(glue("{args$python_source}/unify_groups.py"))
  tryCatch(
    expr = {
      unmatched_only <- tb %>%
        filter(nchar(ProteinGroupId) == 1) %>%
        mutate(Group = "U")
      has_others <- tb %>%
        filter(nchar(ProteinGroupId) > 1) %>%
        as.data.frame() %>%
        py$findNewGroups() %>%
        as_tibble()
      bind_rows(has_others, unmatched_only)
    },
    error = \(cnd) {
      print(reticulate::py_last_error())
      stop("Caught reticulate error")
    }
  )
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
  if (TEST) {
    combined <- combined %>% slice(., 1:2000)
  }

  combined[combined == ""] <- NA
  combined <- correct_ids(combined)
  combined <- combined %>% mutate(
    num_peps = sapply(peptideIds, function(x) {
      return(str_count(x, ";") + 1)
    }, USE.NAMES = FALSE),
    mass = sapply(mass, function(x) {
      if (x == "-" | is.na(x)) {
        return(NA)
      } else {
        return(as.double(x))
      }
    }, USE.NAMES = FALSE)
  )

  combined <- combined %>%
    mutate_all(~ replace(., . == "-", NA)) %>%
    mutate_all(~ replace(., . == NaN, NA)) %>%
    mutate(
      GO_evidence = apply(., 1, getEvidence),
      length = as.double(gsub("unknown", NA, length)),
      unique_peptides = lapply(peptideIds, \(x) {
        x <- str_split_1(x, ";") %>%
          lapply(., clean_peptide) %>%
          unlist()
        return(paste0(unique(x), collapse = ";"))
      }) %>% unlist(),
      num_unique_peps = sapply(unique_peptides, \(x) {
        return(str_count(x, ";") + 1)
      }, USE.NAMES = FALSE),
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
  if (MAP_PFAMS && "GO" %in% colnames(combined)) {
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
        combined$GO_max_sv <- gc$GO_max_sv
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
  if (MAP_KEGG && any(grepl("KEGG", colnames(combined)))) {
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

  ## Merge with quantification data
  if (MERGE_QUANT) {
    cat("BEGIN: merging quantification\n", file = LOGFILE, append = TRUE)
    combined <- merge_with_quant(combined, read_tsv(args$directlfq), "directlfq") %>%
      mutate(across(contains("directlfq"), base::log2))
    combined <- merge_with_quant(combined, read_tsv(args$maxlfq), "maxlfq")
    combined <- get_top3(args$directlfq_aqreformat, combined)
    cat("COMPLETE: merging quantification\n", file = LOGFILE, append = TRUE)
  }

  ## Find new groups
  if (UNIFY) {
    cat("BEGIN: unifying groups\n", file = LOGFILE, append = TRUE)
    combined <- unify_groups(combined)
    cat("COMPLETE: unifying groups\n", file = LOGFILE, append = TRUE)
  }

  ## Arrange columns
  rename_lookup <- c(
    eggNOG_preferred_name = "Preferred_name",
    eggNOG_description = "Description"
  )
  combined <- get_organism(combined, denovo_org = str_replace_all(args$denovo_org, "_", " "))
  combined <- combined %>%
    dplyr::rename(any_of(rename_lookup)) %>%
    relocate(where(is.numeric),
      .after = where(is.character)
    ) %>%
    relocate(c("q.value", "posterior_error_prob"), .before = "q_adjust") %>%
    relocate(any_of(c("peptideIds", "seq")), .after = where(is.numeric)) %>%
    relocate(contains("interpro")) %>%
    relocate(contains("KEGG")) %>%
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
  parser <- add_option(parser, "--directlfq", type = "character")
  parser <- add_option(parser, "--directlfq_aqreformat", type = "character")
  parser <- add_option(parser, "--maxlfq", type = "character")
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
