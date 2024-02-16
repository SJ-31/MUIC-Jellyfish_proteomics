library(tidyverse)
library(AnnotationForge)

prepDf <- function(annotation_file, path_to_extra) {
  # Format tibble for use with makeOrgPackage function
  # - Split GO column
  # - Add GO evidence column
  # - Merge with proteins downloaded from UniProt intended for comparison
  tb <- read_tsv(annotation_file) %>%
    filter(!is.na(GO)) %>%
    dplyr::select(ProteinId, GO_IDs, GO_evidence, header) %>%
    dplyr::rename(GID = ProteinId)
  compare <- list.files(path_to_extra, "*_reviewed.tsv",
                        full.names = TRUE
  ) %>%
    read_tsv() %>%
    dplyr::select(Entry, GO_IDs, `Protein names`) %>%
    dplyr::filter(!is.na(GO_IDs)) %>%
    dplyr::rename(GID = Entry, header = `Protein names`)
  tb <- bind_rows(tb, compare)
  go_df <- tb %>%
    apply(., 1, function(x) {
      gos <- str_split_1(x["GO_IDs"], "\\||;|,") %>%
        lapply(., gsub, pattern = "_.*", replacement = "") %>%
        unlist()
      evidence <- x["GO_evidence"]
      if (!is.na(evidence)) {
        evidence <- str_split_1(evidence, ",") %>%
          lapply(., gsub, pattern = ":.*", replacement = "") %>%
          unlist()
      } else {
        evidence <- "IEA"
      }
      if (length(evidence) < length(gos)) {
        evidence <- c(evidence, rep("IEA", length(gos) - length(evidence)))
      }
      stacked <- as.list(x) %>%
        lapply(., function(x) {
          rep(x, length(gos))
        }) %>%
        as_tibble() %>%
        mutate(GO = gos, EVIDENCE = evidence)
      return(stacked)
    }) %>%
    bind_rows() %>%
    dplyr::select(c("GID", "GO", "EVIDENCE")) %>%
    as.data.frame() %>%
    unique()
  info_df <- data.frame(
    GID = tb$GID, HEADER = tb$header, ENTREZID = tb$GID
  ) %>% unique()
  return(list(info = info_df, go = go_df))
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-c", "--combined_annotations"))
  parser <- add_option(parser, c("-p", "--path_to_extra"))
  parser <- add_option(parser, c("-a", "--author"))
  parser <- add_option(parser, c("-m", "--maintainer"))
  parser <- add_option(parser, c("-t", "--tax_id"))
  parser <- add_option(parser, c("-s", "--species"))
  parser <- add_option(parser, c("-g", "--genus"))
  args <- parse_args(parser)

  dfs <- prepDf(args$combined_annotations, args$path_to_extra)
  makeOrgPackage(
    gene_info = dfs$info, go = dfs$go,
    version = "0.1",
    maintainer = args$author,
    outputDir = ".",
    author = args$maintainer,
    tax_id = args$tax_id,
    genus = args$genus,
    species = args$species,
    goTable = "go",
    verbose = TRUE
  )
}
