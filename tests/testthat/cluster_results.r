library(tidyverse)
library(fpc)

res <- read_tsv("./cluster_stats.tsv")
leiden <- read_tsv("./leiden_metrics.tsv")
