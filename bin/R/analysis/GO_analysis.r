library(glue)
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))
# Fill this with stuff from go_tests.r


if (sys.nframe() == 0 && length(commandArgs(TRUE))) {
  library("optparse")
  parser <- OptionParser()
  parser <- add_option(parser, c("-r", "--r_source"))
  parser <- add_option(parser, c("-r", "--python_source"))
  parser <- add_option(parser, c("-c", "--combined_results"))
  parser <- add_option(parser, c("-g", "--ontologizer_path"))
  parser <- add_option(parser, c("-o", "--output_path"))
  args <- parse_args(parser)

  d <- goData(args$combined_results,
              onto_path = args$ontologizer_path)
  # Gene Set enrichment analysis
  # But how to resolve ties?
  # Could you adjust by the number of peptides or coverage?
  gene_sets <- list(
    unknown_to_db = d$sample_tb %>%
      filter(inferred_by == "interpro" |
               inferred_by == "eggNOG" |
               grepl("[UDT]", ProteinId)) %>%
      pluck("ProteinId"),
    has_mods = d$sample_tb %>%
      filter(ID_method == "open" | !is.na(Mods)) %>%
      pluck("ProteinId"),
    toxins = names(getToxinProteins(d$prot_go_map$sample))
  )
  pwy <- groupPathways(d$sample_tb)
  gene_sets <- c(gene_sets, pwy$grouped)

  by_intensity <- mergeLfq(d$sample_tb, "mean")
  f <- fgseaWrapper("log_intensity", distinct(by_intensity), gene_sets)
  fgsea_intensity <- plotFgsea(gene_sets, f$ranked, f$result)
}


