options(
  browser = "firefox",
  rlang_backtrace_on_error = "full",
  error = rlang::entrace
)
rlang::global_entrace()
library("reticulate")
library("tidyverse")
library("glue")

main <- function(args) {
  source(glue("{args$r_source}/helpers.r"))
  ont <- new.env()
  reticulate::source_python(glue("{args$python_source}/ontologizer_wrapper.py"), envir = ont)
  combined <- read_tsv(args$input)

  O <- reticulate_show_error(ont$Ontologizer(combined, args$executable, args$go_path))
  groups <- list()
  groups[["id_with_open"]] <- dplyr::filter(combined, ID_method == "open" |
    ID_method == "both") |> pluck("ProteinId")
  # Modified proteins or identified in open search
  groups[["unknown_to_db"]] <- dplyr::filter(combined, inferred_by == "interpro" |
    inferred_by == "eggNOG" |
    grepl("[UDT]", ProteinId)) |> pluck("ProteinId")
  # Proteins not known to database, inferred with eggNOG and interpro
  groups[["unmatched_only"]] <- dplyr::filter(combined, grepl("U", ProteinId) | Group == "U") |> pluck("ProteinId")
  params <- list(`-m` = "Bonferroni-Holm")
  print(groups)
  # Proteins mapped ONLY by unmatched peptides or are themselves unmatched peptides
  results <- O$runAll(groups, params) |> lapply(as_tibble)
  return(results)
}

getSlims <- function(args) {
  # Must also report the number of terms that couldn't be slimmed
  source(glue("{args$r_source}/GO_helpers.r"))
  source(glue("{args$r_source}/helpers.r"))
  ontologizer <- get_ontologizer(args$results_path)
  ont_vectors <- ontologizer %>% discard(\(x) any(str_detect("data.frame", class(x))))
  ont_slims <- sapply(names(ont_vectors), \(.) NULL)
  ont_vectors <- ont_vectors %>% lapply(., \(x) x[1:1000])
  for (group in names(ont_vectors)) {
    all_slims <- lapply(names(ont_vectors[[group]]), \(x) {
      slims <- get_go_slim(x, args$go_path, args$go_slim_path)
      return(slims)
    }) %>% unlist()
    ont_slims[[group]] <- all_slims %>%
      unique() %>%
      discard(is.na)
    unknown_to_db <- go_info_tb(ont_slims$unknown_to_db_GO)
    id_with_open <- go_info_tb(ont_slims$id_with_open_GO)
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
    info_tb <- go_info_tb(go_vec)
    tb <- tb %>%
      mutate(sorted_p = transformP(`p.adjusted`)) %>%
      inner_join(., info_tb, by = join_by(x$ID == y$GO_IDs))
    tb
  }
  source(glue("{args$r_source}/GO_helpers.r"))
  source(glue("{args$r_source}/GO_text_mining_helpers.r"))
  source(glue("{args$r_source}/helpers.r"))
  results <- get_ontologizer(args$results_path)
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
  parser <- add_option(parser, c("-p", "--python_source"))
  parser <- add_option(parser, c("--results_path"))
  parser <- add_option(parser, c("--executable"))
  parser <- add_option(parser, c("-s", "--go_slim_path"))
  parser <- add_option(parser, c("-g", "--go_path"))
  parser <- add_option(parser, c("-m", "--mode"))
  parser <- add_option(parser, c("-d", "--go_tm_dir"))
  args <- parse_args(parser)
  if (args$mode == "prep") {
    m <- main(args)
    lmap(m, \(x) {
      write_tsv(x[[1]], glue("ontologizer-{names(x)}.txt"))
    })
  } else if (args$mode == "get_slims") {
    getSlims(args)
  } else if (args$mode == "word_cloud") {
    wordClouds(args)
  }
}
