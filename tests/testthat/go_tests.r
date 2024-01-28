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
comparison_taxa_files <- list.files("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa",
                                    "*_reviewed.tsv", full.names = TRUE)
comparison_tibs <- comparison_taxa_files %>%
  lapply(\(x) {
    read_tsv(x) %>%
      mutate(`Gene Ontology IDs` = unlist(lapply(
        `Gene Ontology IDs`, gsub, pattern = " ", replacement = "", .))
      ) %>%
      dplyr::rename(GO_IDs = `Gene Ontology IDs`)
  }) %>%
  `names<-`(lapply(comparison_taxa_files, gsub,
                   pattern = ".*/(.*)\\.tsv", replacement = "\\1"))

go_tib_comparison <- names(comparison_tibs) %>%
  lapply(., \(x) {
    tib <- comparison_tibs[[x]] %>% filter(!is.na(GO_IDs))
    gos <- goVector(df = tib, go_column = "GO_IDs")
    return(tibble(GO_IDs = gos, taxon = gsub("_reviewed", "", x)))
  }) %>%
  bind_rows()


ontologizer <- ontoResults(onto_path)
cur_df <- readr::read_tsv(test_path) %>%
  filter(!is.na(GO)) %>%
  mutate(GO_IDs = cleanGO(GO))


sample_name <- "C_indra"
sample_gos <- goVector(cur_df, go_column = "GO_IDs") %>% unique()
all_gos <- c(go_tib_comparison$GO_IDs, sample_gos) %>% unique()
interpro_gos <- goVector(cur_df, go_column = "GO_IDs", "Anno_method", "interpro")
eggnog_gos <- goVector(cur_df, go_column = "GO_IDs", "Anno_method", "eggNOG")
go_tib_comparison <- bind_rows(go_tib_comparison, tibble(GO_IDs = sample_gos,
                                                         taxon = sample_name))
go_tib_sample <- tibble(GO_IDs = sample_gos,
                        sig_from_downloaded_db = sample_gos %in% ontologizer$gos_from_downloads,
                        sig_id_with_open = sample_gos %in% ontologizer$id_with_open)

## Visualizing semantic similarity
# Can calculate GO similarity between proteins, either as go lists
# GOSemSim::mgoSim(toy_list, toy_list2[1:5], semData = semdata$MF)
# or proteins themselves

# Protein clusters
# GOSemSim::geneSim(one[1], one[2], semData = semdata$MF)
# GOSemSim::clusterSim(one, two, semData = semdata$MF)

distances <- GOSemSim::mgeneSim(cur_df$ProteinId,
                                semData = semdata$MF, combine = "BMA",
                                drop = "NULL", measure = "Lin"
)
# Obtain pairwise similarity between all proteins in list

rrvgo <- reduceGOList(sample_gos)
anno_method_pcoa <- pcoaWithTib(distances = distances,
                                cur_df, "ProteinId") %>%
  biplotCustom(., colour_column = "Anno_method", x = "PC1", y = "PC2")
id_method_pcoa <- pcoaWithTib(distances = distances, cur_df,
                              "ProteinId") %>%
  biplotCustom(., colour_column = "ID_method", x = "PC1", y = "PC2")
go_pcoa <- sapply(ONTOLOGIES, \(x) NULL)
for (ontology in ONTOLOGIES) {
  go_pcoa[[ontology]] <- sapply(c("from_downloaded_db", "id_with_open"), \(x) NULL)
  for (interest in names(go_pcoa[[ontology]])) {
    go_pcoa[[ontology]][[ interest ]] <- pcoaWithTib(rrvgo$sim_matrix[[ontology]],
                                                     go_tib_sample, "GO_IDs")
  }
}

# Working with go embeddings
# Two sets of embeddings are being used here: "sample", which contains
# GO terms/proteins found from the proteins identified in the shotgun
# pipeline (i.e. from the C. indra data) and "all", which also
# contain GO terms/proteins of other taxa used for comparison


# Obtain embeddings from python
py$wanted_gos <- all_gos
py_run_string("wanted_gos = set(wanted_gos)")
all_embeddings <- py$loadEmbeddings(embeddings_path, "embds", py$wanted_gos)
sample_embeddings <- all_embeddings %>% keep(names(.) %in% sample_gos)
embedding_tib_all <- all_embeddings %>%
  as_tibble() %>%
  t() %>%
  m2Tib(., "GO_IDs")


prot_tib_sample <- dplyr::select(cur_df, c("ProteinId", "GO_IDs"))
prot_tib_all <- comparison_tibs %>%
  bind_rows() %>%
  dplyr::select(Entry, GO_IDs) %>%
  rename(Entry = "ProteinId")

embedding_tib_sample <- embedding_tib_all %>% filter(GO_IDs %in% sample_gos)

embedding_pca <- embedding_tib_all %>%
  dplyr::select(-GO_IDs) %>%
  rda()

sample_umap <- list(data = umapAndJoin(embedding_tib_sample, go_tib_sample))
sample_umap$plot$from_downloaded_db <- biplotCustom(sample_umap$data,
                                                    "sig_from_downloaded_db")
sample_umap$plot$id_with_open <- biplotCustom(sample_umap$data,
                                              "sig_id_with_open")
comparison_umap <- umapAndJoin(embedding_tib_all, go_tib_comparison)

embeddings_only_all <- dplyr::select(embedding_tib_all, -GO_IDs)
embeddings_only_sample <- dplyr::select(embedding_tib_sample, -GO_IDs)
# Distance matrices
# "euclideanDistance" and "cosineSimilarity" defined in a2v.py
sample_dm_euclidean <- dist(t2Df(embedding_tib_sample, "GO_IDs"),
                            method = "euclidean")
all_dm_euclidean <- dist(t2Df(embedding_tib_all, "GO_IDs"),
                         method = "euclidean")
sample_dm_cosine <- makeDistMatrix(sample_embeddings,
                                   sample_embeddings, cosineSimilarity)
all_dm_cosine <- makeDistMatrix(all_embeddings, all_embeddings, cosineSimilarity)

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
