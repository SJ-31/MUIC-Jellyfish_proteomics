library(tidyverse)
library(AnnotationForge)

## Just use the custom ProteinIds column instead
## getGeneId <- function(row) {
##   if (!is.na(row["UniProtKB_ID"])) {
##     return(paste0("UniProtKB:", row["UniProtKB_ID"]))
##   }
##   if (!is.na(row["NCBI_ID"])) {
##     return(paste0("NCBI:", row["NCBI_ID"]))
##   }
##   return(paste0("Custom:", row["ProteinId"]))
## }

prepDf <- function(annotation_file) {
  df <- read_tsv(annotation_file)
  has_gos <- df %>%
    filter(!is.na(GO)) %>%
    mutate(GID = ProteinId)
  go_df <- has_gos %>%
    apply(., 1, function(x) {
      gos <- str_split_1(x["GO"], "\\||;") %>%
        lapply(., gsub, pattern = "_.*", replacement = "") %>%
        unlist()
      evidence <- str_split_1(x["GO_evidence"], ",") %>%
        lapply(., gsub, pattern = ":.*", replacement = "") %>%
        unlist()
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
    GID = has_gos$GID, HEADER = has_gos$header, ENTREZID = has_gos$GID
  ) %>% unique()
  return(list(info = info_df, go = go_df))
}

if (sys.nframe() == 0) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-c", "--combined_annotations"))
  parser <- add_option(parser, c("-a", "--author"))
  parser <- add_option(parser, c("-m", "--maintainer"))
  parser <- add_option(parser, c("-t", "--tax_id"))
  parser <- add_option(parser, c("-s", "--species"))
  parser <- add_option(parser, c("-g", "--genus"))
  args <- parse_args(parser)

  dfs <- prepDf(args$combined_annotations)
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
