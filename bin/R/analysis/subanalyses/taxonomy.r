library("waffle")
library("treemap")
library("ggrepel")
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}

PALETTE <- "ggthemes::Classic_20"
TABLES <- list()
GRAPHS <- list()

data <- M$data |>
  inner_join(M$taxa_tb, by = join_by(ProteinId))
toxins <- data |> filter(if_any(contains("GO_category"), \(x) x == "toxin_related"))

cnidaria <- data |> filter(Phylum == "Cnidaria")
others <- data |> filter(Phylum != "Cnidaria")

GRAPHS$phyla_coverage <- gg_numeric_dist(list(
  cnidaria$pcoverage_align,
  others$pcoverage_align
), method = "boxplot") + xlab("Phylum") + ylab("Coverage (%)") +
  guides(fill = element_blank()) + scale_x_discrete(labels = c("Cnidaria", "Other"))


# Check if Cnidarian proteins have higher coverage (they should)
tests <- wilcox.test(cnidaria$pcoverage_align, others$pcoverage_align, alternative = "greater") |>
  to("data.name", "Coverage of Cnidarian proteins against other taxa") |>
  to("alternative", "Cnidarian proteins have greater coverage") |>
  htest2tb()
tests <- bind_rows(
  tests,
  wilcox.test(cnidaria$num_unique_peps, others$num_unique_peps,
    alternative = "greater"
  ) |>
    to("data.name", "Number of unique peptides identified from Cnidarian proteins against other taxa") |>
    to("alternative", "Cnidarian proteins have more unique peptides") |>
    htest2tb()
)

cubozoa <- cnidaria |> filter(Class == "Cubozoa")

# Sanity check for origin of proteins from other taxa
infer_table <- with(data, {
  table(Phylum == "Cnidaria", inferred_by == "initial_database")
})
chi <- infer_table |> chisq.test()
or <- get_odds_ratio(infer_table)
# The odds of proteins from other taxa being inferred by the initial database search
# is `or` times Cnidarian proteins.
# We expect that it should be less, since the database consists primarily of Cnidarian
# proteins

# Check roles of non-Cnidarian proteins
# TODO: Update this when you have COG categories
# list(others$GO_category_MF, cnidaria$GO_category_MF) |>
#   lapply(\(x) {
#     table <- lapply(x, str_split_1, pattern = ";") |>
#       unlist() |>
#       table() |>
#       table2tb(id_col = "GO_category_MF")
#   })

# What are the proteins in other taxa that have high coverage?
highest_coverage <- others |>
  arrange(desc(pcoverage_align))
highest_coverage$header

# ----------------------------------------
# Relationship between taxonomic status and matching
tb <- M$data |>
  select(ProteinId, header, MatchedPeptideIds) |>
  mutate(
    matched_by_T = grepl("T", MatchedPeptideIds),
    matched_by_D = grepl("D", MatchedPeptideIds),
    matched_by_U = grepl("U", MatchedPeptideIds)
  ) |>
  inner_join(M$taxa_tb, by = join_by(ProteinId)) |>
  filter(!is.na(Class))


matched_by <- c("matched_by_T", "matched_by_D", "matched_by_U")
classes <- tb$Class |> unique()

class_tests <- chisqNME(tb, matched_by, "Class", "matched_by", "Class")

TABLES$class_tests_chi <- class_tests$gt$chi
TABLES$class_tests_contingency <- class_tests$gt$class_tests_contingency




save(c(TABLES, GRAPHS), glue("{M$outdir}/taxonomy"))
