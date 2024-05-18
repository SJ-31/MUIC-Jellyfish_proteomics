library(tidyverse)
library(glue)
library(fpc)

if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  setwd("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/testthat")
} else {
  setwd("/mnt/data/shannc/nf/tests/testthat")
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}


clusterResults <- function(path) {
  return(
    list(
      members = read_tsv(glue("{path}/cluster_members.tsv")),
      stats = read_tsv(glue("{path}/cluster_stats.tsv")),
      leiden = read_tsv(glue("{path}/leiden_metrics.tsv"))
    )
  )
}


sortHelper <- function(metric, tb, ascending = TRUE) {
  selected <- select(tb, all_of(c(metric, "method", "parameter")))
  if (ascending) {
    return(arrange(selected, !!as.symbol(metric)))
  }
  return(arrange(selected, desc(!!as.symbol(metric))))
}


findBest <- function(stat_tb) {
  metrics <- list()
  # Metrics to maximize
  metrics$dunn <- sortHelper("dunn", stat_tb, ascending = FALSE)
  metrics$calinhara <- sortHelper("calinhara", stat_tb, ascending = FALSE)
  #
  # Average silhouette score
  metrics$avg_sil <- sortHelper("silhouette_width", stat_tb, ascending = FALSE)

  # Correlation between distances and a 0-1-vector where 0 means same cluster,
  # 1 means different clusters

  # Metrics to minimize
  # The average similarity measure of each cluster with its most similar
  # cluster
  metrics$davies_bouldin <- sortHelper("davies_bouldin", stat_tb, ascending = TRUE)

  # The average distance within clusters

  # The widest distance within clusters

  # Ratio of average within-cluster distances / average between-cluster distances

  winners <- lapply(names(metrics), \(x) {
    current <- metrics[[x]]
    winner <- slice(current, 1)
    return(glue("{winner$method} {winner$parameter}"))
  }) %>% `names<-`(names(metrics))
  return(list(
    winners = winners,
    metrics = metrics
  ))
}

source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin/"),
  figure_path = glue("{wd}/tests/testthat/output/figures"),
  combined_results = glue("{wd}/results/C_indra_A/1-First_pass/C_indra_all.tsv"),
  embedding_path = glue("{wd}/tests/nf-test-out/C_indra_a2v_go_embeddings/distances.hdf5"),
  sample_name = "C_indra"
)

dpath <- "../nf-test-out/C_indra_prottrans_embeddings/distances.hdf5"
epath <- "../nf-test-out/C_indra_prottrans_embeddings/embeddings.hdf5"
e <- embeddingData(
  args$combined_results,
  args$sample_name,
  epath,
  dpath
)

OUTDIR <- glue("{wd}/tests/testthat/output/")


a2v <- clusterResults(glue("{OUTDIR}/cluster_go_a2v"))
sem <- clusterResults(glue("{OUTDIR}/cluster_go_semantic"))
pt <- clusterResults(glue("{OUTDIR}/cluster_prottrans"))


a2v_best <- findBest(a2v$stats)
sem_best <- findBest(sem$stats)
pt_best <- findBest(pt$stats)
results <- purrr::reduce(
  list(a2v_best$winners, sem_best$winners, pt_best$winners),
  \(x, y) bind_rows(as_tibble(x), as_tibble(y))
) %>% mutate(embedding_source = c("a2v", "semantic", "prottrans"))
sem_best
a2v_best


# 2024-05-11
# the best performing are variations of hclust. But
# hclust with "average" linkage is the top scorer for 7/12 tests
# so we'll use that. But note that the best cut height is rather inconsistent
# between metrics
# Best cut heights:
# Prottrans: 1/10 of the maximum height of 0.1
# Wang semantic sim: mostly inconsistent between metrics
# a2v: (1/10 doing well is likely an outlier, nothing joined up). Winning metrics for this cannot be trusted. Do not use a2v for clustering then
