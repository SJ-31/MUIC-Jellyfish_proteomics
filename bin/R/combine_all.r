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
  "ProteinId", "ProteinGroupId", "header", "category",
  "Group", "matchedPeptideIds", "NCBI_ID", "UniProtKB_ID",
  "organism", "lineage", "ID_method", "Anno_method"
)

HEADER_QUERIES <- list(
  unknown = "unknown|uncharacterized|predicted|unnamed|hypothetical",
  venom_component = "toxin|porin",
  transport = "pump|transport|channel",
  membrane = "membrane",
  translation = "ribosom|elongation|initiation|translation",
  cytoskeleton = "actin|tubulin|cytoskele|dynein|kinesin|catenin",
  other = ""
)


#' Assign a protein into a loose category
#'
#' @description
#' Protein is categorized according to its header description
headerCategorize <- function(header) {
  lvec <- map_lgl(names(HEADER_QUERIES), \(x) {
    grepl(HEADER_QUERIES[[x]], header, ignore.case = TRUE)
  })
  return(names(HEADER_QUERIES)[lvec][1])
}

##'   Remove description from a vector of GO terms separated by ";"
##'   tODO: change the splitting pattern to ; once you have data afterwards
cleanGO <- function(go_vector) {
  go_vector %>%
    lapply(., \(x) {
      if (is.na(x)) {
        return(x)
      }
      s <- str_split_1(x, pattern = ",|;") %>% keep(grepl("GO:", .))
      r <- unlist(lapply(s, gsub, pattern = "_.*", replacement = "")) %>% unique()
      return(paste0(r, collapse = ";"))
    }) %>%
    unlist()
}

##'   Splits an entry with multiple matched peptides so that a match
##'   attempt to the quantification tibble is made with every peptide
splitPeptidesMatch <- function(tb, quant_tb) {
  split_up <- tibbleDuplicateAt(tb, "matchedPeptideIds", ";") %>%
    nest(.by = ProteinId) %>%
    apply(., 1, \(x) {
      group_tb <- x[["data"]]
      id <- x[["ProteinId"]]
      duplicated <- group_tb[1,] %>% mutate(matchedPeptideIds = id)
      return(bind_rows(duplicated, group_tb) %>%
               mutate(
                 ProteinId = id,
                 .before = "ProteinGroupId"
               ))
    }) %>%
    bind_rows()
  return(left_join(split_up, quant_tb, by = join_by(x$matchedPeptideIds == y$ProteinId)))
}


##' When there are multiple quantified peptides for a protein, take
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
      joined_ids <- paste0(x[["data"]]$matchedPeptideIds, collapse = ";")
      id <- x[["ProteinId"]]
      return(top_three[1,] %>%
               select(-contains(quant_name)) %>%
               mutate(.,
                      matchedPeptideIds = joined_ids,
                      ProteinId = id
               ) %>%
               bind_cols(means))
    }) %>%
    bind_rows()
  return(tb)
}


mergeWithQuant <- function(main_tb, quant_tb, quant_name) {
  # Merge annotation tibble with quantification tibble
  # Entries that have been matched by multiple peptides
  # peptides (de novo, transcriptome etc.) in blast and interpro are
  # extracted and handled differently:
  # Compute average and median values across samples
  has_multiple <- main_tb %>% filter(!is.na(matchedPeptideIds) &
                                       ProteinId != matchedPeptideIds)
  if (nrow(has_multiple) != 0) {
    has_multiple <- splitPeptidesMatch(has_multiple, quant_tb)
    has_multiple <- meanTop3(has_multiple, quant_name = quant_name)
  }
  full_proteins <- main_tb %>% filter(is.na(matchedPeptideIds) |
                                        ProteinId == matchedPeptideIds)
  full_proteins <- left_join(full_proteins, quant_tb,
                             by = join_by(x$ProteinId == y$ProteinId)
  )
  bound <- bind_rows(full_proteins, has_multiple)
  calcs <- bound %>%
    dplyr::rowwise() %>%
    dplyr::reframe(
      "{quant_name}_mean" := mean(c_across(contains(quant_name))),
      "{quant_name}_median" := median(c_across(contains(quant_name)))
    )
  return(bind_cols(bound, calcs))
}

writeAlignments <- function(row, file_name) {
  header <- ifelse(row[["header"]] == "unknown", row[["ProteinId"]],
                   row[["header"]]
  )
  pep_count <- row[["num_unique_peps"]]
  write.fasta(row[["seq"]], header, open = "a", file.out = file_name)
  write.fasta(row[["alignment"]],
              glue("ALIGNED PEPTIDES | COUNT: {pep_count}"),
              open = "a",
              file.out = file_name
  )
  cat("\n", file = file_name, append = TRUE)
}

organismFromHeader <- function(row) {
  # Parse the organism from an NCBI or uniprot ID
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
  "lineage", "GO", "GO_evidence", "KEGG_Genes", "PANTHER", "matchedPeptideIds",
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
  v <- str_split_1(values, ";") %>% as.double()
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
                                  #' Parse GO evidence codes for all entries that have them
                                  #' GOs added by eggNOG and interpro are automatically labelled IEA
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
      if (x == "-" | is.na(x)) {
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
    select(-names(redundant[redundant])) %>%
    mutate_all(~replace(., . == "-", NA)) %>%
    mutate_all(~replace(., . == NaN, NA)) %>%
    mutate(
      GO_evidence = apply(., 1, getEvidence),
      length = as.double(gsub("unknown", NA, length)),
      unique_peptides = lapply(peptideIds, \(x) {
        x <- str_split_1(x, ";") %>%
          lapply(., cleanPeptide) %>%
          unlist()
        return(paste0(unique(x), collapse = ";"))
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
      pep_adjust = unlist(lapply(posterior_error_prob, sortVals)),
      GO_IDs = cleanGO(GO),
      category = map_chr(header, headerCategorize)
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

  ## Map pfam domains to GO
  source_python(glue("{args$r_source}/interpro_api.py"))
  combined <- py$mapPfams(
    to_annotate = as.data.frame(combined),
    p2g_path = args$pfam2go,
    i2g_path = args$interpro2go,
    pfam_db_path = args$pfam_db
  ) %>% as_tibble()

  ## Merge with quantification data
  directlfq <- read_tsv(args$directlfq)
  flashlfq <- read_tsv(args$flashlfq)
  combined <- mergeWithQuant(combined, directlfq, "directlfq")
  combined <- mergeWithQuant(combined, flashlfq, "flashlfq")

  ## Categorize

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
  parser <- add_option(parser, "--output", type = "character")
  parser <- add_option(parser, "--interpro", type = "character")
  parser <- add_option(parser, "--eggnog", type = "character")
  parser <- add_option(parser, "--downloads", type = "character")
  parser <- add_option(parser, "--fdr", type = "double")
  parser <- add_option(parser, "--alignment_file", type = "character")
  parser <- add_option(parser, "--pep_thresh", type = "double")
  parser <- add_option(parser, "--is_denovo", type = "character")
  parser <- add_option(parser, "--directlfq", type = "character")
  parser <- add_option(parser, "--flashlfq", type = "character")
  parser <- add_option(parser, "--interpro2go", type = "character")
  parser <- add_option(parser, "--pfam2go", type = "character")
  parser <- add_option(parser, "--pfam_db", type = "character")
  parser <- add_option(parser, "--coverage",
                       type = "character",
                       default = TRUE,
                       action = "store_true"
  )
  parser <- add_option(parser, "--empai",
                       type = "character",
                       default = TRUE,
                       action = "store_true"
  )
  parser <- add_option(parser, "--sort_mods",
                       type = "character",
                       default = TRUE,
                       action = "store_true"
  )
  parser <- add_option(parser, "--r_source", type = "character")
  args <- parse_args(parser)
  source(glue("{args$r_source}/helpers.r"))
  results <- main(args)
  write_tsv(results$all, args$output)
}
