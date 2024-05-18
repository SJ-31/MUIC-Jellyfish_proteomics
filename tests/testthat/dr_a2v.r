EMBEDDING_TYPE <- "GO"
if (grepl("Bio_SDD", getwd())) {
  setwd("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/testthat")
  DATA <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data"
} else {
  setwd("/mnt/data/shannc/nf/tests/testthat")
  DATA <- "/mnt/data/shannc/nf/data"
}
source("./test_dr_data.r")


e$metadata <- d$go_tb$sample %>%
  mutate(.,
    sub_ontology =
      map_chr(GO_IDs, termOntology)
  ) %>%
  filter(GO_IDs %in% rownames(e$cosine))
e$color <- "sub_ontology"
join_on <- "GO_IDs"


for (t in techniques) {
  a2v_path <- glue("{args$figure_path}/a2v_all")
  alabel <- labelGen(glue("{EMBEDDING_TYPE}-level {t}"), "test", "Using a2v distances")
  a2v_dr <- drWrapper(e, join_on = join_on, glue("{a2v_path}/{t}"), "a2v", t)
  plotDr(a2v_dr$to_plot, e$color, a2v_path, t, alabel, TRUE)
}
