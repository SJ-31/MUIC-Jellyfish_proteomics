# File for testing out GO analysis

## Load samples
source("./bin/R/GO_helpers.r")
source("./bin/R/rrvgo_modified.r")
sample <- "./results/C_indra/1-First_pass/C_indra_all.tsv"
onto_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/ontologizer/"
orgdb_pth <- "./tests/testthat/output/org.Cindrasaksajiae.eg.db"
rrvgo_path <- "./tests/testthat/output/rrvgo"
db_name <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)
orgdb <- prepOrgDb(orgdb_pth)

semdata <- lapply(ONTOLOGIES, \(x) {
  GOSemSim::godata(OrgDb = db_name, ont = x, keytype = "GID")
}) %>% `names<-`(ONTOLOGIES)
test_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/1-First_pass/jellyfish_all.tsv"
embeddings_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/go_embedded.npz"

sample_name <- "C_indra"

d <- goDataGlobal(
  uniprot_data_dir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa",
  sample_data = test_path,
  sample_name = sample_name,
  onto_path = onto_path
)


#' Visualizing semantic similarity
#'  Can calculate GO similarity between proteins, either as go lists
#' GOSemSim::mgoSim(toy_list, toy_list2[1:5], semData = semdata$MF)
#' or proteins themselves

# Protein clusters
# GOSemSim::geneSim(one[1], one[2], semData = semdata$MF)
# GOSemSim::clusterSim(one, two, semData = semdata$MF)

# Obtain pairwise similarity between all proteins in list
prot_dist_sample <- GOSemSim::mgeneSim(
  genes = dplyr::mutate(d$sample_tb,
                        ProteinId = map_chr(ProteinId, \(x) gsub("-SAMPLE", "", x))
  )$ProteinId,
  semData = semdata$MF, combine = "BMA",
  drop = "NULL", measure = "Wang"
)

go_freq <- d$sample_tb$GO_IDs %>%
  lapply(., str_split_1, pattern = ";") %>%
  unlist() %>%
  table() %>%
  sort(decreasing = TRUE) %>%
  as_tibble() %>%
  rename(c("GO_ID" = ".", "count" = "n")) %>%
  mutate(
    term = map_chr(GO_ID, \(x)
      ifelse(is.null(GOTERM[[x]]), NA, GOTERM[[x]]@Term)),
    ontology = map_chr(GO_ID, \(x)
      ifelse(is.null(GOTERM[[x]]), NA, GOTERM[[x]]@Ontology))
  )

rrvgo_sample <- reduceGOList(d$go_vec$sample)
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

## TODO: Protein set profiling
# groups <- groupGO(
#   gene = cur_df$ProteinId, OrgDb = db_name,
#   ont = "MF", level = 2, keyType = "GID"
# )
# group_df <- as_tibble(groups@result)
#
# pcoa <- as.data.frame(vegan::wcmdscale(distances, k = 2)) %>%
#   tibble::rownames_to_column(var = "ProteinId") %>%
#   as_tibble()
# pcoa <- inner_join(pcoa, cur_df, by = join_by(x$ProteinId == y$ProteinId))
# x <- "V1"
# y <- "V2"
# ggplot(pcoa, aes(x = .data[[x]], y = .data[[y]]))
#
#

gene_sets <- list(
  unknown_to_db = d$sample_tb %>%
    filter(Anno_method == "interpro" | Anno_method == "eggNOG", grepl("[UDT]", ProteinId)) %>%
    pluck("ProteinId"),
  from_open_search = d$sample_tb %>%
    filter(ID_method == "open") %>%
    pluck("ProteinId"),
  toxins = names(getToxinProteins(d$prot_go_map$sample))
) %>% lapply(., \(x) map_chr(x, \(x) gsub("-SAMPLE", "", x)))

# Ranked gene set enrichment analysis
ranked <- d$sample_tb %>%
  dplyr::select(ProteinId, directlfq_mean) %>%
  dplyr::mutate(ProteinId = map_chr(ProteinId, \(x) gsub("-SAMPLE", "", x))) %>%
  dplyr::filter(!is.na(directlfq_mean)) %>%
  column_to_rownames(var = "ProteinId") %>%
  (\(x) {
    y <- as.vector(x$directlfq_mean)
    names(y) <- rownames(x)
    return(y)
  }) %>%
  sort(., decreasing = TRUE)


fgsea(gene_sets, ranked, scoreType = "pos")
