options(
  browser = "firefox",
  rlang_backtrace_on_error = "full",
  error = rlang::entrace
)
rlang::global_entrace()
library(tidyverse)
library(glue)

goVector <- function(df, column, filter) {
  # Return a flattened vector of GO terms from df
  df <- df %>% filter(!!as.symbol(column) == filter)
  gos <- df$GO %>%
    lapply(., strsplit, split = "\\||;|,") %>%
    unlist() %>%
    discard(is.na) %>%
    discard(!grepl("GO:", .))
}

writeOz <- function(file, df) {
  # Write ontologizer-formatted df to file
  header <- "GoStat IDs Format Version 1.0"
  write_lines(header, file)
  write_tsv(df, file, append = TRUE)
}

ozFormat <- function(tb) {
  # Format df as input for ontologizer program
  # Optional filtering by specific columns and values
  prepped <- tb %>%
    select(c(ProteinId, GO)) %>%
    filter(!is.na(GO)) %>%
    mutate(GO = sapply(GO, function(x) {
      s <- unlist(strsplit(x, "\\||;|,")) %>% discard(!grepl("GO:", .))
      x <- sapply(s, gsub,
        pattern = "_.*", replacement = "",
        USE.NAMES = FALSE
      )
      return(paste0(x, collapse = ",")) # comma delimiter is required for ontologizer formatting
    }, USE.NAMES = FALSE)) %>%
    filter(GO != "NA")
  return(prepped)
}

prep <- function(args) {
  combined <- read_tsv(args$input)
  u <- ozFormat(combined)
  # Protein universe
  io <- ozFormat(dplyr::filter(combined, ID_method == "open" |
    ID_method == "both"))
  # Modified proteins or identified in open search
  sa <- ozFormat(
    combined %>% dplyr::filter(inferred_by == "interpro" |
      inferred_by == "eggNOG" |
      grepl("[UDT]", ProteinId))
  )
  # Proteins not known to database, inferred with eggNOG and interpro
  d <- ozFormat(
    dplyr::filter(combined, grepl("U", ProteinId) | Group == "U")
  )
  # Proteins mapped ONLY by unmatched peptides or are themselves unmatched peptides
  return(list(
    universe = u,
    id_open = io,
    standard_annotation = sa,
    unmatched_only = d
  ))
}

getSlims <- function(args) {
  # Must also report the number of terms that couldn't be slimmed
  source(glue("{args$r_source}/GO_helpers.r"))
  source(glue("{args$r_source}/helpers.r"))
  ontologizer <- ontoResults(args$results_path)
  ont_vectors <- ontologizer %>% discard(\(x) any(str_detect("data.frame", class(x))))
  ont_slims <- sapply(names(ont_vectors), \(.) NULL)
  ont_vectors <- ont_vectors %>% lapply(., \(x) x[1:1000])
  for (group in names(ont_vectors)) {
    all_slims <- lapply(names(ont_vectors[[group]]), \(x) {
      slims <- getGoSlim(x, args$go_path, args$go_slim_path)
      return(slims)
    }) %>% unlist()
    ont_slims[[group]] <- all_slims %>%
      unique() %>%
      discard(is.na)
    unknown_to_db <- goInfoTb(ont_slims$unknown_to_db_GO)
    id_with_open <- goInfoTb(ont_slims$id_with_open_GO)
    write_tsv(unknown_to_db, "unknown_to_db_GO_slims.tsv")
    write_tsv(id_with_open, "id_with_open_GO_slims.tsv")
  }
}

wordClouds <- function(args) {
  transformP <- function(p_vec) {
    lg <- -log(p_vec)
    lg[lg == Inf] <- .Machine$integer.max
    lg
  }
  prep <- function(tb, go_vec) {
    info_tb <- goInfoTb(go_vec)
    tb <- tb %>%
      mutate(sorted_p = transformP(`p.adjusted`)) %>%
      inner_join(., info_tb, by = join_by(x$ID == y$GO_IDs))
    tb
  }
  source(glue("{args$r_source}/GO_helpers.r"))
  source(glue("{args$r_source}/GO_text_mining_helpers.r"))
  source(glue("{args$r_source}/helpers.r"))
  results <- ontoResults(args$results_path)
  params <- list(
    term_col = "name", sort_by = "sorted_p", compound = FALSE,
    color_col = "ontology", shape = "circle", word_size = "sorted_p"
  )
  id_with_open <- prep(results$id_with_open, names(results$id_with_open_GO))
  unknown_to_db <- prep(results$unknown_to_db, names(results$unknown_to_db_GO))
  with_open_tk <- tokenize2Plot(id_with_open, params)
  unknown_tk <- tokenize2Plot(unknown_to_db, params)
  with_open_cloud <- wordcloudCustom(
    with_open_tk$tb, params,
    with_open_tk$abbrevs
  )
  unknown_cloud <- wordcloudCustom(unknown_tk$tb, params, unknown_tk$abbrevs)
  ggsave(with_open_cloud, filename = "id_with_open_wordcloud.png")
  ggsave(unknown_cloud, filename = "unknown_to_db_wordcloud.png")
}


if (sys.nframe() == 0 && length(commandArgs(TRUE))) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-i", "--input"))
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("--results_path"))
  parser <- add_option(parser, c("-s", "--go_slim_path"))
  parser <- add_option(parser, c("-g", "--go_path"))
  parser <- add_option(parser, c("-m", "--mode"))
  parser <- add_option(parser, c("-d", "--go_tm_dir"))
  args <- parse_args(parser)
  if (args$mode == "prep") {
    m <- prep(args)
    writeOz("protein_mappings.ids", m$universe)
    writeLines(m$id_open$ProteinId, "id_with_open.txt")
    writeLines(m$unmatched_only$ProteinId, "unmatched_only.txt")
    writeLines(m$standard_annotation$ProteinId, "unknown_to_db.txt")
    writeLines(m$universe$ProteinId, "universe.txt")
  } else if (args$mode == "get_slims") {
    getSlims(args)
  } else if (args$mode == "word_cloud") {
    print(args)
    wordClouds(args)
  }
}
