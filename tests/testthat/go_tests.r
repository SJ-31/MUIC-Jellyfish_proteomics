library(tidyverse)
library(glue)
# File for testing out GO analysis
if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  reticulate::use_virtualenv("/home/shannc/.venv/default")
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
  protrans_dist = glue("{wd}/tests/nf-test-out/C_indra_prottrans_embeddings/distances.hdf5"),
  protrans_embd = glue("{wd}/tests/nf-test-out/C_indra_prottrans_embeddings/embeddings.hdf5"),
  dist_path = glue("{wd}/tests/nf-test-out/C_indra_esm_embeddings/distances.hdf5"),
  go_path = glue("{wd}/data/reference/go.obo"),
  go_slim_path = glue("{wd}/data/reference/goslim_generic.obo")
)

## Load samples
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/rrvgo_modified.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
orgdb_pth <- "./tests/testthat/output/org.Cindrasaksajiae.eg.db"
DB_NAME <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)
rrvgo_path <- "./tests/testthat/output/rrvgo"

uniprot <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/protein_databases/comparison_taxa/reviewed_all.tsv"

sample_name <- "C_indra"


d <- goData(args$combined_results,
  onto_path = args$ontologizer_path
)

e <- embeddingData(
  args$combined_results,
  args$sample_name,
  args$protrans_embd,
  args$protrans_dist,
)

id_with_open <- d$ontologizer$id_with_open_GO
dist <- e$cosine



nested$KEGG_Pathway <- associatedPathways(e$metadata, nested)
nested %>% filter(!is.na(KEGG_Pathway))


# Maybe plot kegg pathways???

chosen <- "ko00512"
sample_tb <- pathway_tbs[[chosen]]
sample <- pathway_graphs[[chosen]]
sample |>
  activate(edges) |>
  as_tibble() |>
  View()

ggraph(sample, layout = "fr") +
  geom_edge_diagonal(aes(color = subtype_name)) +
  geom_node_point(aes(filter = !type %in% c("map", "compound")))

gg <- ggraph(sample, layout = "kk")
gg + geom_edge_link(
  aes(
    color = subtype_name,
  )
) +
  geom_node_point(
    aes(filter = !type %in% c("map", "compound")),
    fill = gg$data[!gg$data$type %in% c("map", "compound"), ]$bgcolor,
    color = "black",
    shape = 21, size = 4
  ) +
  geom_node_point(
    aes(filter = !type %in% c("map", "gene")),
    fill = gg$data[!gg$data$type %in% c("map", "gene"), ]$bgcolor,
    color = "black",
    shape = 21, size = 6
  ) +
  geom_node_text(
    aes(
      label = name,
    ),
    repel = TRUE,
    bg.colour = "white"
  )
