library("xml2")
library("tidygraph")
library("ggkegg")

#' Node vectors
#'
#' @description
#' Return all nodes from a ggkegg igraph object
getNodes <- function(g) {
  if (is.null(g)) {
    return(NA)
  } else {
    return(
      g %>%
        activate(nodes) %>%
        as_tibble() %>%
        pluck("name") %>%
        flattenJoined(., " ")
    )
  }
}

filterKnown <- function(pathway_tb) {
  pathway_tb |> filter((name %in% GENES | type != "gene") &
    (name %in% KO | type != "ortholog") & (name %in% PATHWAYS | type != "map"))
}

pathwayCompleteness <- function(pathway_tb) {
  wanted_types <- c("gene", "ortholog", "map")
  is_known <- filterKnown(pathway_tb)
  relevant <- pathway_tb |> filter(type %in% wanted_types)
  return(nrow(filter(is_known, type %in% wanted_types)) / nrow(relevant))
}

#' Ranked gene set enrichment analysis
#'
#' @description
#' Rank proteins in tb according to a specified column, in
#' descending fashion
#' @param quant the column to rank proteins on
fgseaWrapper <- function(quant, tb, gene_sets, id_col = "ProteinId") {
  ranked <- tb %>%
    dplyr::select(all_of(c(id_col, quant))) %>%
    dplyr::filter(!is.na(!!quant)) %>%
    column_to_rownames(var = id_col) %>%
    (\(x) {
      y <- as.vector(x[[quant]])
      names(y) <- rownames(x)
      return(y)
    }) %>%
    sort(., decreasing = TRUE)
  # Ranks are descending, so most abundant first
  fgsea <- fgsea::fgsea(pathways = gene_sets, ranked, scoreType = "pos")
  # Recommended to switch to "pos"
  return(list(result = fgsea, ranked = ranked))
}

# See here for more plotting
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
plotFgsea <- function(gene_sets, ranked_list, fgsea_result, p_cutoff = 0.05) {
  fgsea_plots <- list()
  for (n in seq_len(nrow(fgsea_result))) {
    if (fgsea_result[n, ]$padj < p_cutoff) {
      set <- fgsea_result[n, ]$pathway
      padjust <- round(fgsea_result[n, ]$padj, 6)
      fgsea_plots[[set]] <- plotEnrichment(gene_sets[[set]], ranked_list) +
        labs(title = glue("Set: {set}, adjusted p-value = {padjust}"))
    }
  }
  return(fgsea_plots)
}

fgseaGroup <- function(tb, grouping_col, gene_sets) {
  group_map <- hash::hash()
  group_map[tb$ProteinId] <- tb[[grouping_col]]
  grouped_intensity <- tb %>%
    group_by(!!as.symbol(grouping_col)) %>%
    nest() %>%
    mutate(
      mean_intensity = map_dbl(data, \(x) mean(x$log_intensity, na.rm = TRUE)),
      sd_intensity = map_dbl(data, \(x) sd(x$log_intensity, na.rm = TRUE)),
      members = lapply(data, \(x) x$ProteinId),
      size = map_dbl(data, \(x) nrow(x))
    ) %>%
    filter(!is.na(mean_intensity) & !is.na(!!as.symbol(grouping_col)))

  gene_sets_groups <- lapply(gene_sets, \(x) map_unique(x, group_map))
  f <- fgseaWrapper("mean_intensity", grouped_intensity, gene_sets_groups, id_col = grouping_col)
  return(list(fgsea = f, groups = gene_sets_groups))
}

get_pathway_title <- function(pathway) {
  pathway_file <- glue("{CACHE}/{pathway}.xml")
  if (!file.exists(pathway_file)) {
    warning(glue("Pathway file for {pathway} not downloaded"))
    return(NULL)
  }
  xml <- read_xml(pathway_file)
  return(xml_attr(xml, "title"))
}

#' Reformat into proteins x term
#'
#' @description
#' Helper function to format tibble into a binary dataframe where rows are proteins and
#' columns are GO term names. Values indicate whether or not a protein has that GO term
#'
#' @param tb long-form tibble containing ProteinIds, their mapped
#' GO terms and the names of the terms
#' @param terms terms to use for columns
idxterm <- function(tb, terms) {
  make_binary_row <- function(x) {
    binary <- map_dbl(terms, \(v) ifelse(v %in% x["data"]$data$term, 1, 0))
    names(binary) <- terms
    binary <- as.list(binary)
    binary$ProteinId <- x["ProteinId"][[1]]
    as_tibble_row(binary)
  }

  nested <- group_by(tb, ProteinId) |> nest()
  formatted <- apply(nested, 1, make_binary_row) |> bind_rows()
  column_to_rownames(formatted, var = "ProteinId")
}
