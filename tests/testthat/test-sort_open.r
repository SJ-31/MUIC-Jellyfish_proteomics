source("./bin/sort_open_searches.r")
args <- list(
  input_directory = "./results/jellyfish/1-First_pass/Open_search/Percolator/",
  seq_header_file = "./results/jellyfish/Databases/seq-header_mappings.tsv",
  r_source = "~/Bio_SDD/MUIC_senior_project/workflow/bin",
  path = "./results/jellyfish/1-First_pass/Open_search/Percolator"
)
source(glue("{args$r_source}/helpers.r"))
m <- main(args)
