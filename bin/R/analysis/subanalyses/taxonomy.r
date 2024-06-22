# Scratch session 2024-06-10-Monday

TABLES <- list()
GRAPHS <- list()

data <- data |>
  inner_join(taxa_tb, by = join_by(ProteinId))
toxins <- data |> filter(if_any(contains("GO_category"), \(x) x == "toxin_related"))

cnidaria <- data |> filter(Phylum == "Cnidaria")
others <- data |> filter(Phylum != "Cnidaria")

GRAPHS$phyla_coverage <- ggplotNumericDist(list(
  cnidaria$pcoverage_nmatch,
  others$pcoverage_nmatch
), method = "boxplot") + xlab("Phylum") + ylab("Coverage (%)")


# Check if Cnidarian proteins have higher coverage (they should)
tests <- wilcox.test(cnidaria$pcoverage_nmatch, others$pcoverage_nmatch, alternative = "greater") |>
  to("data.name", "Coverage of Cnidarian proteins against other taxa") |>
  to("alternative", "Cnidarian proteins have greater coverage") |>
  htest2Tb()
tests <- bind_rows(
  tests,
  wilcox.test(cnidaria$num_unique_peps, others$num_unique_peps,
    alternative = "greater"
  ) |>
    to("data.name", "Number of unique peptides identified from Cnidarian proteins against other taxa") |>
    to("alternative", "Cnidarian proteins have more unique peptides") |>
    htest2Tb()
)

cnidaria$Class |> table()
cubozoa <- cnidaria |> filter(Class == "Cubozoa")

# Sanity check for origin of proteins from other taxa
infer_table <- with(data, {
  table(Phylum == "Cnidaria", inferred_by == "initial_database")
})
chi <- infer_table |> chisq.test()
or <- oddsRatio(infer_table)
# The odds of proteins from other taxa being inferred by the initial database search
# is `or` times Cnidarian proteins.
# We expect that it should be less, since the database


# Check roles of non-Cnidarian proteins
list(others$GO_category_MF, cnidaria$GO_category_MF) |>
  lapply(\(x) {
    table <- lapply(x, str_split_1, pattern = ";") |>
      unlist() |>
      table() |>
      table2Tb(id_col = "GO_category_MF")
  })

# What are the proteins in other taxa that have high coverage?
highest_coverage <- others |>
  arrange(desc(pcoverage_nmatch)) |>
  slice(1:10)
highest_coverage$header
highest_coverage$organism
