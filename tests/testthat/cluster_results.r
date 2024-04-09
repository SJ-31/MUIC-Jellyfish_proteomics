library(tidyverse)
library(fpc)

res <- read_tsv("./tests/testthat/cluster_stats.tsv")
leiden <- read_tsv("./tests/testthat/leiden_metrics.tsv")

sortHelper <- function(metric, tb, ascending = TRUE) {
  selected <- select(tb, all_of(c(metric, "method", "parameter")))
  if (ascending) {
    return(arrange(selected, !!as.symbol(metric)))
  }
  return(arrange(selected, desc(!!as.symbol(metric))))
}

# Metrics to maximize
dunn <- sortHelper("dunn", res, ascending = FALSE)
dunn2 <- sortHelper("dunn2", res, ascending = FALSE)
# Average silhouette score
avg_sil <- sortHelper("avg.silwidth", res, ascending = FALSE)

# Correlation between distances and a 0-1-vector where 0 means same cluster,
# 1 means different clusters
pearsongamma <- sortHelper("pearsongamma", res, ascending = FALSE)


# Metrics to minimize
# The average distance within clusters
average_within <- sortHelper("average.within", res)
# The widest distance within clusters
widest_within <- sortHelper("widestgap", res)
# Ratio of average within-cluster distances / average between-cluster distances
wb_ratio <- sortHelper("wb.ratio", res)
