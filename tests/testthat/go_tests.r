# File for testing out GO analysis

## Load samples
source("./bin/GO.r")
sample <- "./results/jellyfish//1-First_pass/jellyfish_all.tsv"
onto_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/nf-test-out/ontologizer/"
orgdb_pth <- "./tests/testthat/output/org.Cindrasaksajiae.eg.db"
db_name <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)
orgdb <- prepOrgDb(orgdb_pth)

semdata <- lapply(ONTOLOGIES, \(x) {
  GOSemSim::godata(OrgDb = db_name, ont = x, keytype = "GID")
}) %>% `names<-`(ONTOLOGIES)
test_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/1-First_pass/jellyfish_all.tsv"
embeddings_path <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/go_embedded.npz"

sample_name <- "C_indra"


d <- goDataGlobal(uniprot_data_dir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa",
                  sample_data = test_path,
                  sample_name = sample_name,
                  onto_path = onto_path)


#' Visualizing semantic similarity
#'  Can calculate GO similarity between proteins, either as go lists
#'GOSemSim::mgoSim(toy_list, toy_list2[1:5], semData = semdata$MF)
#'or proteins themselves

#Protein clusters
#GOSemSim::geneSim(one[1], one[2], semData = semdata$MF)
#GOSemSim::clusterSim(one, two, semData = semdata$MF)

# Obtain pairwise similarity between all proteins in list
prot_dist_sample <- GOSemSim::mgeneSim(d$sample_tb$ProteinId,
                                       semData = semdata$MF, combine = "BMA",
                                       drop = "NULL", measure = "Lin"
)


toxins <- sample_tb %>% filter(grepl("[Tt]oxin", header))


rrvgo_sample <- reduceGOList(d$go_vec$sample)
rrvgo_all <- reduceGOList(d$go_vec$all)
anno_method_pcoa <- pcoaWithTb(distances = prot_dist_sample,
                               d$sample_tb, "ProteinId") %>%
  biplotCustom(., colour_column = "Anno_method", x = "PC1", y = "PC2")
id_method_pcoa <- pcoaWithTb(distances = prot_dist_sample, d$sample_tb,
                             "ProteinId") %>%
  biplotCustom(., colour_column = "ID_method", x = "PC1", y = "PC2")
go_pcoa <- sapply(ONTOLOGIES, \(x) NULL)

for (ontology in ONTOLOGIES) {
  go_pcoa[[ontology]] <- sapply(c("from_downloaded_db", "id_with_open"), \(x) NULL)
  for (interest in names(go_pcoa[[ontology]])) {
    go_pcoa[[ontology]][[ interest ]] <- pcoaWithTb(rrvgo_sample$sim_matrix[[ontology]],
                                                    d$go_tb$sample, "GO_IDs")
  }
}

## TODO: Protein set profiling
# groups <- groupGO(
#   gene = cur_df$ProteinId, OrgDb = db_name,
#   ont = "MF", level = 2, keyType = "GID"
# )
# group_df <- as_tibble(groups@result)
#
# ## TODO: Ranked gene set enrichment analysis
# pcoa <- as.data.frame(vegan::wcmdscale(distances, k = 2)) %>%
#   tibble::rownames_to_column(var = "ProteinId") %>%
#   as_tibble()
# pcoa <- inner_join(pcoa, cur_df, by = join_by(x$ProteinId == y$ProteinId))
# x <- "V1"
# y <- "V2"
# ggplot(pcoa, aes(x = .data[[x]], y = .data[[y]]))
#
