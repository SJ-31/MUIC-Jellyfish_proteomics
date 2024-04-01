# Working with go embeddings
# Two sets of embeddings are being used here: "sample", which contains
# GO terms/proteins found from the proteins identified in the shotgun
# pipeline (i.e. from the C. indra data) and "all", which also
# contain GO terms/proteins of other taxa used for comparison

d <- goDataGlobal(
  uniprot_data_dir = args$uniprot_data_dir,
  onto_path = args$ontologizer_path,
  sample_data = args$combined_results,
  sample_name = args$sample_name,
  sample_only = SAMPLE_ONLY
)

combineFun <- ifelse(args$protein_embd_mode == "mean",
                     base::mean, base::sum)

embdFromPy <- function(go_vector) {
  reticulate::source_python(glue::glue("{args$python_source}/a2v.py"))
  py$wanted_gos <- go_vector
  py_run_string("wanted_gos = set(wanted_gos)")
  embd <- py$loadEmbeddings(
    args$embeddings_path, "embds",
    py$wanted_gos
  ) %>%
    as_tibble() %>%
    t() %>%
    m2Tb(., first_col = "GO_IDs")
  return(embd)
}

# Obtain embeddings from python
if (!SAMPLE_ONLY) {
  all_embd_go <- embdFromPy(d$go_vec$all)

  d$go_tb$all <- d$go_tb$all %>%
    dplyr::filter(GO_IDs %in% all_embd_go$GO_IDs)

  uniprot_embd <- read_tsv(args$uniprot_embeddings)
  sample_prot_embd <- goEmbedding2Prot(d$prot_go_map$sample, all_embd_go, combineFun)
  all_embd_prot <- bind_rows(uniprot_embd, sample_prot_embd)
  rm(uniprot_embd)

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
    embd_type = args$protein_embedding_mode,
    data = all_embd_prot,
    tb = d$protein$all,
    color = "taxon"
  )
  sample_embd_go <- all_embd_go %>% filter(GO_IDs %in% d$go_vec$sample)
} else {
  sample_embd_go <- embdFromPy(d$go_vec$sample)
  sample_prot_embd <- goEmbedding2Prot(d$prot_go_map$sample, sample_embd_go, combine_func = combineFun)
}
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
  embd_type = args$protein_embedding_mode,
  data = sample_prot_embd,
  tb = d$sample_tb,
  color = c("inferred_by", "ID_method", "category")
)

