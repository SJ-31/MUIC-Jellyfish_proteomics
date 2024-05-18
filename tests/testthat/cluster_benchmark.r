library(fpc)
library(tidyverse)
library(glue)
library(dbscan)
TEST <- FALSE
# Compare clustering GO terms using semantic distance or cosine distances
# with GO embeddings

DIST_TYPE <- "semantic"
EMBEDDING_TYPE <- "go"
# LOG
# 2024-05-06 Finished with dist_type = "a2v" and embedding_type = "go"
# 2024-05-07 Finished with dist_type = "semantic" and embedding_type = "go"
# Finished with embedding_type = "protein"


if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
  env <- "/home/shannc/Bio_SDD/miniconda3/envs/reticulate"
  setwd("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests/testthat")
} else {
  setwd("/mnt/data/shannc/nf/tests/testthat")
  wd <- "/home/shannc/workflow"
  env <- "/home/shannc/anaconda3/envs/reticulate"
}

args <- list(
  r_source = glue("{wd}/bin/R"),
  python_source = glue("{wd}/bin/"),
  figure_path = glue("{wd}/tests/testthat/output/figures"),
  combined_results = glue("{wd}/results/C_indra_A/1-First_pass/C_indra_all.tsv"),
  embedding_path = glue("{wd}/tests/nf-test-out/C_indra_a2v_go_embeddings/distances.hdf5"),
  sample_name = "C_indra"
)

source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/DR_helpers.r"))
source(glue("{args$r_source}/cluster_helpers.r"))
source(glue("{args$r_source}/analysis/prepare_embeddings.r"))
orgdb_pth <- "./output/org.Cindrasaksajiae.eg.db"
db_name <- gsub(".*\\/", "", orgdb_pth, fixed = FALSE)

d <- goData(args$combined_results)
gos <- d$go_vec$sample

if (DIST_TYPE == "semantic" && EMBEDDING_TYPE == "go") {
  prepOrgDb(orgdb_pth)
  dist <- rrvgo::calculateSimMatrix(gos,
    orgdb = db_name,
    ont = "BP", method = "Wang",
    keytype = "GID"
  ) %>%
    simToDist()
} else if (EMBEDDING_TYPE == "go") {
  dpath <- "../nf-test-out/C_indra_a2v_go_embeddings/distances.hdf5"
  epath <- "../nf-test-out/C_indra_a2v_go_embeddings/embeddings.hdf5"
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    epath,
    dpath
  )
  dist <- e$cosine
} else if (EMBEDDING_TYPE == "protein") {
  dpath <- "../nf-test-out/C_indra_prottrans_embeddings/distances.hdf5"
  epath <- "../nf-test-out/C_indra_prottrans_embeddings/embeddings.hdf5"
  e <- embeddingData(
    args$combined_results,
    args$sample_name,
    epath,
    dpath
  )
  dist <- e$cosine
}

if (EMBEDDING_TYPE == "go") {
  OUTDIR <- glue("./output/cluster_go_{DIST_TYPE}")
  JOIN_COL <- "id"
} else {
  OUTDIR <- glue("./output/cluster_prottrans")
  JOIN_COL <- "ProteinId"
}
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR)
}

if (TEST) {
  sample <- rownames(dist) %>% base::sample(size = 500)
  dist <- filterDistMatrix(dist, sample)
  OUTDIR <- glue("./output/cluster_testzone")
}

# decrease size for testing
CLUSTER_STATS <- NULL
CLUSTER_MEMBERS <- NULL

clusterBenchmark <- function(dist, technique, cluster_stats, cluster_members) {
  if (technique == "hclust") {
    param <- "height"
    heights <- seq(min(dist), max(dist), length.out = 12)[-1]
    heights <- round(heights[-length(heights)], 2)
    cluster_labels <- benchmarker(
      dist,
      \(d, h) hclustSk(d, h, linkage = "average"), as.double(heights), param, technique
    )
  } else if (technique == "hclust_structured") {
    param <- "height"
    heights <- seq(min(dist), max(dist), length.out = 12)[-1]
    heights <- round(heights[-length(heights)], 2)
    cluster_labels <- benchmarker(
      dist,
      \(d, h) hclustSk(d, h, linkage = "average", structured = TRUE),
      as.double(heights), param, technique
    )
  } else if (technique == "hclust_complete") {
    param <- "height"
    heights <- seq(min(dist), max(dist), length.out = 12)[-1]
    heights <- round(heights[-length(heights)], 2)
    cluster_labels <- benchmarker(
      dist,
      \(d, h) hclustSk(d, h, linkage = "complete"), as.double(heights), param, technique
    )
  } else if (technique == "hdbscan") {
    param <- "min_points"
    cluster_labels <- benchmarker(dist, `_hdbscan`, seq(5, 30, by = 5), param, technique)
  } else if (technique == "leiden") {
    param <- "quality_fun"
    # # Mask relationships with euclidean distances
    ## greater than the 3rd quartile for leiden (just set them to the maximum distance away)
    dist[dist > quantile(dist, 0.75)] <- max(dist)
    # # Benchmark leiden
    py$graph <- leidenCreateGraph(dist)
    cluster_labels <- benchmarker(
      dist, `_leiden`,
      c("Modularity", "RBER", "RB", "CPM", "Surprise"),
      "partition_type", "leiden"
    )
    write_tsv(LEIDEN_METRICS, glue("{OUTDIR}/leiden_metrics.tsv"))
  } else if (technique == "protein_groups" && EMBEDDING_TYPE != "go") {
    numberProteinGroups <- function(group_vector) {
      t <- table(group_vector)
      groups <- seq_along(t) %>% `names<-`(names(t))
      return(groups)
    }
    param <- "NONE"
    tb <- read_tsv(args$combined_results) %>%
      filter(ProteinId %in% rownames(dist))
    g <- numberProteinGroups(tb$Group)
    cluster_labels <- map_dbl(tb$Group, \(x) g[x]) %>%
      `names<-`(tb$ProteinId) %>%
      discard(is.na)
    dist <- filterDistMatrix(dist, names(cluster_labels))
    cluster_labels <- cluster_labels[!duplicated(names(cluster_labels))]
    cluster_labels <- list(first = cluster_labels)
  }
  members <- saveClusters(
    cluster_labels, technique, param,
    cluster_members, JOIN_COL
  )
  stats <- getClusterMetrics(dist, cluster_labels, technique, cluster_stats, param)
  return(list(stats = stats, members = members))
}


to_test <- c(
  "hclust",
  "hclust_structured",
  "hclust_complete",
  "hdbscan",
  "leiden"
)
if (EMBEDDING_TYPE == "protein") {
  to_test <- c(to_test, "protein_groups")
}


if (all(diag(dist) != 0)) {
  diag(dist) <- 0
}
diag(dist) <- 0

for (t in to_test) {
  c <- clusterBenchmark(dist, t, CLUSTER_STATS, CLUSTER_MEMBERS)
  dir <- glue("{OUTDIR}/{t}")
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  CLUSTER_STATS <- c$stats
  CLUSTER_MEMBERS <- c$members
  write_tsv(CLUSTER_MEMBERS, glue("{dir}/cluster_members.tsv"))
  write_tsv(CLUSTER_STATS, glue("{dir}/cluster_stats.tsv"))
}

previous_stats <- list.files(OUTDIR, pattern = "stats.tsv", recursive = TRUE, full.names = TRUE)
lapply(previous_stats, read_tsv) %>%
  bind_rows() %>%
  distinct() %>%
  write_tsv(glue("{OUTDIR}/cluster_stats.tsv"))

previous_members <- list.files(OUTDIR,
  pattern = "members.tsv",
  recursive = TRUE, full.names = TRUE
)
seen_cols <- c()
num_rows <- NULL
lapply(previous_members, \(x) {
  tsv <- read_tsv(x)
  if (is.null(num_rows)) {
    num_rows <<- nrow(tsv)
  }
  not_seen <- colnames(tsv)[!colnames(tsv) %in% seen_cols]
  seen_cols <<- c(not_seen, seen_cols)
  if (nrow(tsv) > num_rows) {
    tsv <- dplyr::slice(tsv, 1:num_rows)
  } else if (nrow(tsv) < num_rows) {
    tsv <- tsv %>% add_row()
  }
  return(dplyr::select(tsv, all_of(not_seen)))
}) %>%
  bind_cols() %>%
  write_tsv(glue("{OUTDIR}/cluster_members.tsv"))
