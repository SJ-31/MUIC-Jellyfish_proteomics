library(glue)
# Working with go and protein-langauge model embeddings
# Two sets of embeddings are being used here: "sample", which contains
# GO terms/proteins found from the proteins identified in the shotgun
# pipeline (i.e. from the C. indra data) and "compare", which also
# contain GO terms/proteins of other taxa used for comparison

# Combine GO terms assigned to each protein in protein_map using
# combine_func
goEmbedding2Prot <- function(protein_map, embedding_tb, combine_func) {
  purrr::map(names(protein_map), ~ {
    dplyr::filter(embedding_tb, GO_IDs %in% protein_map[[.x]]) %>%
      reframe(across(where(is.numeric), combine_func)) %>%
      mutate(ProteinId = .x, .before = V1)
  }) %>%
    bind_rows() %>%
    dplyr::filter(!if_any(contains("V"), is.na)) # bug: This shouldn't be necessary
}


a2vGetGOEmbd <- function(go_vector) {
  reticulate::source_python(glue::glue("{args$python_source}/a2v.py"))
  py$wanted_gos <- go_vector
  py_run_string("wanted_gos = set(wanted_gos)")
  embd <- py$loadEmbeddings(
    args$embedding_path, "embds",
    py$wanted_gos
  ) %>%
    as.data.frame() %>%
    t()
  rownames(embd) <- map_chr(
    rownames(embd),
    \(x) gsub("\\.", ":", x)
  )
  colnames(embd) <- paste0("V", seq(dim(embd)[2]))
  return(as.data.frame(embd))
}

getProtEmbd <- function(embd_path, dist_path) {
  reticulate::source_python(glue::glue("{args$python_source}/get_distances.py"))
  results <- getSaved(embd_path, dist_path)
}

combineGOEmbd <- function(embd_path, go_vector, map, write_distances) {
  reticulate::source_python(glue::glue("{args$python_source}/get_distances.py"))
  # Retrieve embeddings per GO id
  embd_go <- a2vGetGOEmbd(go_vector) %>%
    rownames_to_column(., var = "GO_IDs") %>%
    as_tibble()
  combined <- goEmbedding2Prot(map, embd_go, base::mean)
  names <- combined$ProteinId
  embeddings <- combined %>%
    dplyr::select(-ProteinId) %>%
    as.matrix()
  writeEmbeddingsHDF5(embd_path, names, embeddings)
  if (!missing(write_distances) && !is.null(write_distances)) {
    writeDistances(embeddings, names, write_distances)
  }
}

embeddingData <- function(
    combined_results,
    sample_name,
    embedding_path,
    # Path to embeddings in hdf5
    dist_path, # Path to distance matrix in hdf5
    comparison_meta # Add this argument when
    # comparing sample proteins against external
    ) {
  py_embd <- getProtEmbd(embedding_path, dist_path)
  py_embd$euclidean <- as.matrix(py_embd$euclidean)
  py_embd$cosine <- as.matrix(py_embd$cosine)
  data <- read_tsv(combined_results) %>%
    distinct(ProteinId, .keep_all = TRUE) %>%
    filter(ProteinId %in% rownames(py_embd$embeddings))
  if (!missing(comparison_meta)) {
    comp_meta <- read_tsv(comparison_meta) %>% rename(ProteinId = Entry)
    data <- bind_rows(data, comp_meta) %>%
      mutate(Taxon = sample_name) %>%
      select(c(ProteinId, Taxon)) %>%
      filter(ProteinId %in% rownames(py_embd$embeddings))
    color <- "Taxon"
  } else {
    color <- c("inferred_by", "ID_method", "GO_category_CC", "GO_category_MF")
  }
  return(list(
    embd = py_embd$embeddings,
    euclidean = py_embd$euclidean,
    cosine = py_embd$cosine,
    metadata = data,
    color = color
  ))
}


#' When called from the command line, the script extracts
#' anc2vec's pre-downloaded GO embeddings (not protein)
#' If `go_only`, then the embeddings of the GO terms in the sample
#' are retrieved. Otherwise, the embeddings for GO terms are combined
#' for a protein-level representation using the mean
if (sys.nframe() == 0 && length(commandArgs(TRUE))) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, "--r_source", type = "character")
  parser <- add_option(parser, "--python_source", type = "character")
  parser <- add_option(parser, "--sample_tsv", type = "character")
  parser <- add_option(parser, "--embd_output", type = "character", help = "output file when writing embeddings")
  parser <- add_option(parser, "--embedding_path",
    type = "character",
    help = "path to file containing embeddings for the proteins in sample_tsv"
  )
  parser <- add_option(parser, "--dist_output", type = "character")
  parser <- add_option(parser, "--comparison_embd", type = "character")
  parser <- add_option(parser, "--comparison_tsv", type = "character")
  parser <- add_option(parser, "--go_only", action = "store_true", default = FALSE)
  parser <- add_option(parser, "--sample_name", type = "character")

  args <- parse_args(parser)
  source(glue("{args$r_source}/GO_helpers.r"))
  if (args$go_only) {
    d <- get_go_data(args$sample_tsv)
    gos <- a2vGetGOEmbd(d$go_vec$sample)
    reticulate::source_python(glue::glue("{args$python_source}/get_distances.py"))
    names <- rownames(gos)
    writeDistances(gos, names, args$dist_output)
    writeEmbeddingsHDF5(args$embd_output, names, as.matrix(gos))
  } else if (is.null(args$comparison_embd)) {
    d <- get_go_data(args$sample_tsv)
    combineGOEmbd(args$embd_output, d$go_vec$sample, d$prot_go_map$sample, args$dist_output)
  } else {
    # If compare, filter only proteins of interest
    id_col <- "ProteinId"
    reticulate::source_python(glue::glue("{args$python_source}/get_distances.py"))
    d <- get_go_data(args$sample_tsv,
      uniprot_tsv = args$uniprot_tsv,
      sample_name = args$sample_name
    )
    sample_embd <- hdf5ToDf(args$embedding_path)
    comparison_embd <- hdf5ToDf(args$comparison_embd)
    # Retrieve embeddings

    # Apply filtering criteria
    sample_meta <- d$sample_tb %>% filter(category == "venom_component")
    min <- min(sample_meta$length)
    comparison_meta <- d$protein$compare %>%
      filter(Length >= min) %>%
      adjustProteinNumbers(nrow(sample_meta), .)

    combined <- comparison_embd %>%
      df2tb(., id_col) %>%
      filter(ProteinId %in% comparison_meta$ProteinId) %>%
      bind_rows(., df2tb(sample_embd, id_col))

    names <- combined$ProteinId
    embd <- column_to_rownames(., id_col)
    writeDistances(names, embd, args$dist_output)
    writeEmbeddingsHDF5(args$embd_output, names, as.matrix(embd))
  }
}
