if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
# --------------------------------------------------------
# Investigating trends in missing quantification
current <- M$data
cq <- c("directlfq", "maxlfq", "flashlfq")
missing_quant_tests <- list()
for (q in cq) {
  noq <- current %>% filter(is.na(!!as.symbol(glue("{cq}_mean"))))
  hasq <- current %>% filter(!ProteinId %in% noq$ProteinId)
  missing_quant_tests[[q]] <- wilcox.test(noq$num_peps, hasq$num_peps,
    alternative = "l"
  )
}
capture.output(missing_quant_tests, file = glue("{M$outdir}/missing_quantification_tests.txt"))
rm(noq)
rm(hasq)
# --------------------------------------------------------

# Show that property of every protein in a Percolator group being matched
# by the exact same peptides gets lost when creating new groups via union-find
data <- M$data
nested <- data |>
  group_by(Group) |>
  nest() |>
  mutate(
    unique_peptides = map_dbl(data, \(x) {
      x$peptideIds |>
        unique() |>
        length()
    }),
    size = map_dbl(data, \(x) nrow(x))
  ) |>
  arrange(desc(size))

# -------------------------------------------------------------

# Confirm that cases where the peptides are longer than the proteins they were
# matched to are the situation where proteins are fragmentary
isFragment <- Logical() ? function(header = ? Character()) {
  str_detect(header, "([fF]ragment)|(partial)")
}
has_longer <- data |> filter(maxPeptideLength(peptideIds) > length)
others <- data |> filter(!ProteinId %in% has_longer)
percent <- sum(map_lgl(has_longer$header, isFragment)) / sum(map_lgl(others$header, isFragment)) * 100
fragments <- data |> filter(map_lgl(header, isFragment))


# ----------------------------------------
# Identify enriched terms by intensity
filterIntensity <- function(tb, quantile = ? Character()) {
  if (quantile == "first") {
    filterFun <- \(x) filter(x, x$log_intensity <= quantile(x$log_intensity, 0.25))
  } else if (quantile == "second") {
    filterFun <- \(x) {
      filter(x, x$log_intensity > quantile(x$log_intensity, 0.25) &
        x$log_intensity < quantile(x$log_intensity, 0.75))
    }
  } else if (quantile == "third") {
    filterFun <- \(x) filter(x, x$log_intensity >= quantile(x$log_intensity, 0.75))
  }
  tb |>
    filterFun() |>
    purrr::pluck("ProteinId")
}

if (!file.exists(glue("{M$ontologizer_path}/high_intensity.tsv"))) {
  ont <- new.env()
  reticulate::source_python(glue("{M$python_source}/ontologizer_wrapper.py"), envir = ont)
  by_intensity <- mergeLfq(data, "mean") %>%
    inner_join(., dplyr::select(data, ProteinId, GO_IDs), by = join_by(ProteinId)) |>
    filter(!is.na(log_intensity))
  O <- ont$Ontologizer(by_intensity, M$ontologizer_exec, M$go_path)
  groups <- list(
    low_intensity = filterIntensity(by_intensity, "first"),
    medium_intensity = filterIntensity(by_intensity, "second"),
    high_intensity = filterIntensity(by_intensity, "third")
  )
  params <- list(`-m` = "Bonferroni-Holm")
  enriched_intensity <- O$runAll(groups, params)
  enriched_intensity |> names()
  lmap(enriched_intensity, \(x) {
    write_tsv(x[[1]], glue("{M$ontologizer_path}/{names(x)}.tsv"))
  })
  plot <- ggplot(by_intensity, aes(x = log_intensity)) +
    geom_histogram(fill = "#69d2e7") +
    xlab("log 10 intensity")
  ggsave(glue("{M$outdir}/figures/intensity_histogram.svg"), plot)
} else {
  intensity <- lapply(c("low", "medium", "high"), \(x) read_tsv(glue("{M$ontologizer_path}/{x}_intensity.tsv")))
  # TODO analyze differences between these
}
