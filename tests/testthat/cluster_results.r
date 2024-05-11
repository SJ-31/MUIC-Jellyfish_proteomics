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

OUTDIR <- glue("{wd}/tests/testthat/output/")

a2v <- clusterResults(glue("{OUTDIR}/cluster_go_a2v"))
sem <- clusterResults(glue("{OUTDIR}/cluster_go_semantic"))
pt <- clusterResults(glue("{OUTDIR}/cluster_prottrans"))


a2v_best <- findBest(a2v$stats)
sem_best <- findBest(sem$stats)
pt_best <- findBest(pt$stats)
# 2024-05-05
# Looks like hierarchichal clustering at 1/10 of the maximum distance
# performs the best
# But this doesn't really make sense...
# For the a2v "clusters", there is literally only one cluster with
# three members. Could this be some form of outlier?
