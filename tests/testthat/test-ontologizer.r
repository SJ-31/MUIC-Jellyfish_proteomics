source("./bin/ontologizer.r")
args <- list(
  input = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/1-First_pass/C_indra_all_wcoverage.tsv",
  python_source = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin",
  executable = "/home/shannc/Bio_SDD/tools/Ontologizer.jar",
  go_path = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/go.obo"
)
m <- main(args)
