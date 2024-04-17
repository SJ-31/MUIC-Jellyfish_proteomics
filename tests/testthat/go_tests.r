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
  uniprot_data_dir = glue("{wd}//data/protein_databases/comparison_taxa"),
  combined_results = glue("{wd}//results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv"),
  ontologizer_path = glue("{wd}//tests/nf-test-out/ontologizer/"),
  embedding_path = glue("{wd}/data/reference/go_embedded.npz"),
  dist_path = glue("{wd}//tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"
  ))

## Load samples
source("./bin/R/GO_helpers.r")
source("./bin/R/analysis/metric_functions.r")
source("./bin/R/rrvgo_modified.r")
sample <- "./results/C_indra/1-First_pass/C_indra_all.tsv"
orgdb_pth <- "./tests/testthat/output/org.Cindrasaksajiae.eg.db"
db_name <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)
rrvgo_path <- "./tests/testthat/output/rrvgo"
# orgdb <- prepOrgDb(orgdb_pth)
# semdata <- lapply(ONTOLOGIES, \(x) {
#   GOSemSim::godata(OrgDb = db_name, ont = x, keytype = "GID")
# }) %>% `names<-`(ONTOLOGIES)

uniprot <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa/reviewed_all.tsv"

sample_name <- "C_indra"


d <- goData(args$combined_results,
            onto_path = args$ontologizer_path)
sample_tb <- read_tsv(args$combined_results)


# Protein clusters
# GOSemSim::geneSim(one[1], one[2], semData = semdata$MF)
# GOSemSim::clusterSim(one, two, semData = semdata$MF)

# go_freq <- d$sample_tb$GO_IDs %>%
#   lapply(., str_split_1, pattern = ";") %>%
#   unlist() %>%
#   table() %>%
#   sort(decreasing = TRUE) %>%
#   as_tibble() %>%
#   rename(c("GO_ID" = ".", "count" = "n")) %>%
#   mutate(
#     term = map_chr(GO_ID, \(x)
#       ifelse(is.null(GOTERM[[x]]), NA, GOTERM[[x]]@Term)),
#     ontology = map_chr(GO_ID, \(x)
#       ifelse(is.null(GOTERM[[x]]), NA, GOTERM[[x]]@Ontology))
#   )
#


# rrvgo_sample <- reduceGOList(d$go_vec$sample)
# # rrvgo_sample$reduced_matrix %>% lmap()
# rrvgo_all <- reduceGOList(d$go_vec$all)
# anno_method_pcoa <- pcoaWithTb(
#   distances = prot_dist_sample,
#   d$sample_tb, "ProteinId"
# ) %>%
#   biplotCustom(., colour_column = "Anno_method", x = "PC1", y = "PC2")
# id_method_pcoa <- pcoaWithTb(
#   distances = prot_dist_sample, d$sample_tb,
#   "ProteinId"
# ) %>%
#   biplotCustom(., colour_column = "ID_method", x = "PC1", y = "PC2")
# go_pcoa <- sapply(ONTOLOGIES, \(x) NULL)
#
# for (ontology in ONTOLOGIES) {
#   go_pcoa[[ontology]] <- sapply(c("from_downloaded_db", "id_with_open"), \(x) NULL)
#   for (interest in names(go_pcoa[[ontology]])) {
#     go_pcoa[[ontology]][[interest]] <- pcoaWithTb(
#       rrvgo_sample$sim_matrix[[ontology]],
#       d$go_tb$sample, "GO_IDs"
#     )
#   }
# }

flattenJoined <- function(vec, split) {
  vec %>%
    discard(is.na) %>%
    lapply(., \(x) str_split_1(x, split)) %>%
    unlist()
}

## TODO: Protein set profiling
groups <- groupGO(
  gene = d$sample_tb$ProteinId, OrgDb = db_name,
  ont = "MF", level = 2, keyType = "GID"
)

sample <- read_tsv("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/combined/UNSPECIFIED_all.tsv")
tb <- sample
library(ggkegg)
spwy <- "ko01502"
ALL_KOS <- flattenJoined(tb$KEGG_ko, "[;,]") %>%
  unique() %>%
  paste0("ko:", .)
ALL_PATHWAYS <- flattenJoined(tb$KEGG_Pathway, "[;,]") %>%
  unique() %>%
  paste0("path:", .)
g <- ggkegg::pathway(spwy, directory = ".temp")
nodes <- g %>%
  activate(nodes) %>%
  as_tibble() %>%
  pluck("name") %>%
  flattenJoined(., " ")
percent_found <- sum(nodes %in% c(ALL_KOS, ALL_PATHWAYS)) / length(nodes)
# This percentage isn't necessarily indicative of the pathway completeness,
# as some pathway components could be optional or swapped out for others

# Filter graph so that only nodes that were found get retained
# Then
filtered <- g %>%
  activate(nodes) %>%
  filter(map_lgl(name, \(x) {
    split <- str_split_1(x, " ")
    return(any(split %in% c(ALL_KOS, ALL_PATHWAYS)))
  }))

gg <- ggraph(filtered, layout = "manual", x = x, y = y) +
  geom_node_rect(fill = "red", aes(filter = type == "ortholog")) +
  overlay_raw_map(spwy) +
  theme_void()


# Check module completeness
ALL_MODULES <- flattenJoined(tb$KEGG_Module, "[;,]") %>% unique()
check_complete <- lapply(ALL_MODULES, \(x) {
  mk <- tryCatch(
    expr = ggkegg::module(x),
    error = \(cnd) NULL
  )
  if (!is.null(mk)) {
    completeness <- ggkegg::module_completeness(mk, ALL_KOS)
    completeness$module <- mk@ID
    completeness$name <- mk@name
    return(completeness)
  } else {
    return(NULL)
  }
})

