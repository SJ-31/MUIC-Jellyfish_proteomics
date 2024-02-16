library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)
args <- list(r_source = "./bin")
source(glue("{args$r_source}/helpers.r"))
final_results <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/1-First_pass/jellyfish_all.tsv"
percolator_dir <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/jellyfish/1-First_pass/Percolator"
percolator_files <- list.files(percolator_dir,
  pattern = "_percolator_proteins.tsv", full.names = TRUE
)
final_results <- read_tsv(final_results)

engines <- percolator_files %>%
  map_chr(., \(x) gsub("_percolator_proteins.tsv", "", x) %>%
    gsub(".*/", "", .))
names(percolator_files) <- engines
percolator_files <- as.list(percolator_files)
tibbles <- percolator_files %>%
  lmap(., \(x) {
    tibbleDuplicateAt(read_tsv(x[[1]]), "ProteinId", ",") %>% mutate(engine = names(x))
  })

id_list <- lapply(engines, function(x) {
  return(filter(bind_rows(tibbles), engine == x)$ProteinId)
}) %>% `names<-`(engines)
svg("engine_intersections.svg")
venn(id_list,
  plotsize = 100, borders = FALSE, box = FALSE,
  ilcs = 1, sncs = 1.5,
  zcolor = "red, blue, green, white, orange, yellow"
)

filtered_tibbles <- lmap(tibbles, \(x) inner_join(x[[1]],
  final_results,
  by = join_by("ProteinId")
))

# TODO: Extract different columns (e.g. lineage, GO) from "filtered tibbles" to see
# if there are statistically significant differences between the different engines
# Are the engines biased in how they find proteins??
