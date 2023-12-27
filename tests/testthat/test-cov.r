library(testthat)
library(glue)
bin <- "../../bin"
source(glue("{bin}/protein_coverage.r"))

intersect <- "../toy_results/toy_intersect.tsv"
map <- "../toy_results/all_normal_mapping.tsv"
prot_df <- read.delim(intersect, sep = "\t") %>% as_tibble() %>% select(-c(use.names, coverage, sequence))
mapping <- read.delim(map, sep = "\t") %>% as_tibble()
mapping <- setNames(as.list(mapping$seq), mapping$id)
covered <- coverage_calc(prot_df, mapping)
write_delim(covered, "~/coverage.sfs", delim = "\t")
