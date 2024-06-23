#' Are there any meaningful differences between the sub-groups of proteins
#' identified in the data?
#'

e <- embeddingData(
  M$combined_results,
  M$sample_name,
  M$prottrans_embd,
  M$prottrans_dist
)


data <- e$metadata %>% select(., -matches("flashlfq|directlfq|maxlfq"))
#' Try to predict what an unannotated protein does based on what the members
#' of its group does (e.g. its GO_IDs and KEGG)

data %>% colnames()

data$ID_method %>% table()
