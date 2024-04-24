library(tidyverse)
library(glue)
# File for testing out GO analysis
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
} else {
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin"),
  embd_type = "protein",
  sample_name = "C_indra",
  uniprot_data_dir = glue("{wd}/data/protein_databases/comparison_taxa"),
  combined_results = glue("{wd}/results/C_indra_A/1-First_pass/C_indra_all_wcoverage.tsv"),
  ontologizer_path = glue("{wd}/tests/nf-test-out/ontologizer/"),
  embedding_path = glue("{wd}/data/reference/go_embedded.npz"),
  dist_path = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"
  ))

## Load samples
source("./bin/R/GO_helpers.r")
source("./bin/R/analysis/metric_functions.r")
library(ggkegg)
library(KEGGREST)

flattenJoined <- function(vec, split) {
  vec %>%
    discard(is.na) %>%
    lapply(., \(x) str_split_1(x, split)) %>%
    unlist()
}

# Goal: determine which pathways are most complete and worth visualizing
ALL_KOS <- flattenJoined(tb$KEGG_ko, "[;,]") %>%
  unique()
ALL_PATHWAYS <- flattenJoined(tb$KEGG_Pathway, "[;,]") %>%
  unique()

# For testing,
ALL_PATHWAYS <- ALL_PATHWAYS[1:10]

pathway_graphs <- vector(mode = "list", length = length(ALL_PATHWAYS)) %>%
  `names<-`(ALL_PATHWAYS)
for (pwy in names(pathway_graphs)) {
  graph <- NA
  try(graph <- ggkegg::pathway(pwy, directory = ".temp"))
  if (!all(is.na(graph))) {
    pathway_graphs[[pwy]] <- graph
  }
}
ALL_KOS_prefixed <- addPrefix(ALL_KOS, "ko:")
ALL_PATHWAYS_prefixed <- addPrefix(ALL_PATHWAYS, "path:")

pathway_nodes <- lapply(pathway_graphs, getNodes)

# Completeness if the ratio of pathways found in the sample against the number
# of nodes in the graph object
# This percentage isn't necessarily indicative of the pathway completeness,
# as some pathway components could be optional or swapped out for others
pathway_completeness <- tibble(
  pathway = ALL_PATHWAYS,
  completeness = map_dbl(pathway_nodes, \(x) {
    if (all(is.na(x))) return(NA)
    return((sum(x %in% c(ALL_KOS_prefixed, ALL_PATHWAYS_prefixed)) / length(nodes)) * 100)
  })
)


# Filter graph so that only nodes that were found get retained
# Then
filtered <- g %>%
  activate(nodes) %>%
  filter(map_lgl(name, \(x) {
    split <- str_split_1(x, " ")
    return(any(split %in% c(paste0("ko:", ALL_KOS), paste0("pathway:", ALL_PATHWAYS))))
  }))

gg <- ggraph(filtered, layout = "manual", x = x, y = y) +
  geom_node_rect(fill = "red", aes(filter = type == "ortholog")) +
  overlay_raw_map(spwy) +
  theme_void()


# Check module completeness
ALL_MODULES <- flattenJoined(tb$KEGG_Module, "[;,]") %>% unique()
ALL_MODULES <- lapply(ALL_MODULES, \(x) {
  return(tryCatch(
    expr = ggkegg::module(x, use_cache = TRUE),
    error = \(cnd) NULL
  ))
}) %>% `names<-`(ALL_MODULES)

check_complete <- ALL_MODULES %>%
  lmap(., \(x) {
    mk <- x[[1]]
    if (is.null(mk)) return(list(NULL))
    completeness <- ggkegg::module_completeness(mk, ALL_KOS)
    completeness$module <- mk@ID
    completeness$name <- mk@name
    return(completeness)
  }) %>%
  discard(is.null) %>%
  bind_rows() %>%
  filter(present_num > 0)

found <- ALL_MODULES %>% discard(is.null)
found$M00001 %>%
  module_text() %>%
  plot_module_text() +
  geom_node_rect()


# Get pathway metadata for summarizing the proteome i.e. discussing wha
# t
# pathway categories there are but also to assess pathway completeness
# Completeness could be based on which modules are present?
keggGetSlowly <- slowly(\(x) keggGet(x), rate = rate_delay(0.5))
# Essential or the API will block you
pathway_info <- lapply(ALL_PATHWAYS, \(x) {
  get <- NULL
  Sys.sleep(0.5)
  try(expr = get <- keggGetSlowly(x))
  return(get[[1]])
}) %>% discard(is.null)
names(pathway_info) <- lapply(pathway_info, \(x) x)


# You can use this to get a bunch of entries, if the entry fails it won't appear
# in the list
test <- keggGet("hsa01100")
