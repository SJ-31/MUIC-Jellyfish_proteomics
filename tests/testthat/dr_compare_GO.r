EMBEDDING_TYPE <- "go"
if (str_detect(getwd(), "Bio_SDD")) {
  setwd("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/testthat")
  DATA <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data"
} else {
  setwd("/mnt/data/shannc/nf/tests/testthat")
  DATA <- "/mnt/data/shannc/nf/data"
}

source("./test_dr_data.r")

SEMDATA <- lapply(ONTOLOGIES, \(x) {
  GOSemSim::godata(OrgDb = db_name, ont = x, keytype = "GID", computeIC = FALSE)
}) %>% `names<-`(ONTOLOGIES)

if (EMBEDDING_TYPE == "protein") {
  e$color <- "category"
} else {
  go_meta <- goInfoTb(d$go_vec$sample, go_file, go_slim_file) %>%
    filter(!is.na(ontology))
  e$color <- "ontology"
}


# Is limited to only one ontology at a time...
if (EMBEDDING_TYPE == "protein") {
  sem_dist_matrix <- GOSemSim::mgeneSim(d$sample_tb$ProteinId,
    SEMDATA$MF,
    measure = "Wang", combine = "BMA"
  ) %>%
    simToDist()
} else {
  sem_dist_matrix <- rrvgo::calculateSimMatrix(go_meta$id,
    orgdb = db_name,
    ont = "MF",
    method = "Wang",
    keytype = "GID"
  ) %>%
    simToDist()
}

if (EMBEDDING_TYPE == "go") {
  join_on <- "id"
} else {
  join_on <- "ProteinId"
}

d$sample_tb <- d$sample_tb %>% filter(ProteinId %in% rownames(sem_dist_matrix))
e$embd <- e$embd[rownames(e$embd) %in% rownames(sem_dist_matrix), ]
e$metadata <- e$metadata %>% filter(ProteinId %in% rownames(sem_dist_matrix))

e$euclidean <- filterDistMatrix(e$euclidean, rownames(sem_dist_matrix))
e$cosine <- filterDistMatrix(e$cosine, rownames(sem_dist_matrix))

sem <- list(
  dist = as.data.frame(sem_dist_matrix),
  metadata = d$sample_tb,
  color = e$color
)
if (EMBEDDING_TYPE == "go") {
  e$metadata <- go_meta %>% filter(id %in% rownames(sem$dist))
  sem$metadata <- e$metadata
}

PLOT <- TRUE
if (PLOT) {
  for (t in techniques) {
    sem_path <- glue("{args$figure_path}/semantic")
    a2v_path <- glue("{args$figure_path}/a2v")
    slabel <- labelGen(glue("{EMBEDDING_TYPE}-level {t}"), "test", "Using semantic distances")
    alabel <- labelGen(glue("{EMBEDDING_TYPE}-level {t}"), "test", "Using a2v distances")

    sem_dr <- drWrapper(sem, join_on = join_on, glue("{sem_path}/{t}"), "sem", t)
    plotDr(sem_dr$to_plot, e$color, sem_path, t, slabel, TRUE)

    a2v_dr <- drWrapper(e, join_on = join_on, glue("{a2v_path}/{t}"), "a2v", t)
    plotDr(a2v_dr$to_plot, e$color, a2v_path, t, alabel, TRUE)
  }
}
