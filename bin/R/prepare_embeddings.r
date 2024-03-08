# Working with go embeddings
# Two sets of embeddings are being used here: "sample", which contains
# GO terms/proteins found from the proteins identified in the shotgun
# pipeline (i.e. from the C. indra data) and "all", which also
# contain GO terms/proteins of other taxa used for comparison

d <- goDataGlobal(
  uniprot_data_dir = args$uniprot_data_dir,
  onto_path = args$ontologizer_path,
  sample_data = args$combined_results,
  sample_name = args$sample_name
)

# Obtain embeddings from python
py$wanted_gos <- d$go_vec$all
py_run_string("wanted_gos = set(wanted_gos)")
all_embd_go <- py$loadEmbeddings(
  embeddings_path, "embds",
  py$wanted_gos
) %>%
  as_tibble() %>%
  t() %>%
  m2Tb(., first_col = "GO_IDs")
d$go_tb$all <- d$go_tb$all %>% dplyr::filter(GO_IDs %in% all_embd_go$GO_IDs)

uniprot_embd <- read_tsv(args$uniprot_embd)
sample_prot_embd <- goEmbedding2Prot(d$prot_go_map$sample, all_embd_go, mean)
all_embd_prot <- bind_rows(uniprot_embd, sample_prot_embd)
rm(uniprot_embd)

sample_embd_go <- all_embd_go %>% filter(GO_IDs %in% d$go_vec$sample)
sample_go <- list(
  data = sample_embd_go,
  tb = dplyr::filter(
    d$go_tb$sample,
    is.element(
      GO_IDs,
      sample_embd_go$GO_IDs
    )
  ),
  color = c("sig_downloaded_db", "sig_id_w_open")
)
sample_protein <- list(
  embd_type = protein_embd_mode,
  data = sample_prot_embd,
  tb = d$sample_tb,
  color = c("Anno_method", "ID_method", "category")
)
all_go <- list(
  data = all_embd_go,
  tb = dplyr::filter(
    d$go_tb$all,
    is.element(
      GO_IDs,
      all_embd_go$GO_IDs
    )
  ),
  color = "taxon"
)
all_protein <- list(
  embd_type = protein_embd_mode,
  data = all_embd_prot,
  tb = d$protein$all,
  color = "taxon"
)
