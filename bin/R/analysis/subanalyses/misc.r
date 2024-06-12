# --------------------------------------------------------
# Investigating trends in missing quantification
current <- run$first
cq <- c("directlfq", "maxlfq", "flashlfq")
missing_quant_tests <- list()
for (q in cq) {
  noq <- current %>% filter(is.na(!!as.symbol(glue("{cq}_mean"))))
  hasq <- current %>% filter(!ProteinId %in% noq$ProteinId)
  missing_quant_tests[[q]] <- wilcox.test(noq$num_peps, hasq$num_peps,
    alternative = "l"
  )
}
capture.output(missing_quant_tests, file = glue("{OUTDIR}/missing_quantification_tests.txt"))
rm(noq)
rm(hasq)
# --------------------------------------------------------

# Show that property of every protein in a Percolator group being matched
# by the exact same peptides gets lost when creating new groups via union-find
data <- run$first
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
