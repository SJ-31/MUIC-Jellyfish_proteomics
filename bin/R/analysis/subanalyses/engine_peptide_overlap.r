if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
PALETTE <- "ggthemes::colorblind"
PALETTE2 <- "ggthemes::Classic_10_Medium"
library("ggplot2")
library("ggVennDiagram")
TABLES <- list()
GRAPHS <- list()

open_search_engines <- c("metamorpheusGPTMD", "msfraggerGPTMD", "msfraggerGlyco")
percolator_all <- read_tsv(M$percolator_all)
ENGINES <- percolator_all$engine |> unique()
alignment_types <- c("denovo", "transcriptome", "database", "unmatched_peptide")
standard_search_engines <- ENGINES[!ENGINES %in% open_search_engines]

data <- read_tsv(M$data_path)

num_peptides_matched <- ta$get_engine_counts(M$percolator_all, data) |>
  as_tibble() |>
  distinct()

# #' Group engines that identify the same peptide groups using Jaccard distance

enginesXPeptides <- num_peptides_matched %>%
  tb_transpose()

engine_dist <- vegan::vegdist(enginesXPeptides, method = "jaccard")
dist_longer <- engine_dist |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column(var = "engine") |>
  pivot_longer(cols = -engine)
pcoa <- vegan::wcmdscale(engine_dist, eig = TRUE)$points %>%
  as.data.frame() |>
  rownames_to_column(var = "engine")
GRAPHS$engine_sim_jaccard <- dist_longer |> ggplot(aes(x = engine, y = name, fill = value)) +
  geom_tile() +
  ylab("Engine") +
  xlab("Engine") +
  scale_fill_paletteer_c("ggthemes::Orange-Gold", name = "Jaccard distance") +
  theme(axis.text.x = element_text(angle = 90)) +
  M$default_theme
GRAPHS$engine_sim_jaccard_biplot <- pcoa |> ggplot(aes(x = Dim1, y = Dim2, color = engine)) +
  geom_point()

standard_engine_tbs <- lapply(standard_search_engines, \(x) {
  percolator_all |> filter(engine == x)
}) %>%
  `names<-`(standard_search_engines)

#' Venn diagram for overlap
id_list <- lapply(standard_engine_tbs, \(x) x$ProteinId) %>%
  `names<-`(names(standard_engine_tbs))
venn <- ggVennDiagram(id_list, label = "none") +
  scale_fill_paletteer_c("ggthemes::Classic Blue")
GRAPHS$engine_venn <- venn

save(c(GRAPHS, TABLES), glue("{M$outdir}/engine_characteristics"))
