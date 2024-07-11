library("plotly")
library("waffle")
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
list(others$GO_category_MF, cnidaria$GO_category_MF) |>
  lapply(\(x) {
    table <- lapply(x, str_split_1, pattern = ";") |>
      unlist() |>
      table() |>
      table2tb(id_col = "GO_category_MF")
  })

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


# ----------------------------------------
# Toxin protein analysis


# Use for header and eggnog description
toxin_words <- c(
  "porin", "hemolysin", "phospholipase", "enterotoxin", "disintegrin", "nucleotidase", "secapin", "defensin", "venom",
  "reprolysin", "crisp", "metalloprotease", "cftx", "metalloproteinase", "pest_crys",
  "metallopeptidase", "conotoxin", "serine protease", "metalloendopeptidase",
  "neurotoxin", "latrotoxin"
)

descriptor_cols <- c(
  "PFAMs", "header", "PANTHER", "eggNOG_description",
  "interpro_description"
)

toxin_query <- paste0(toxin_words, collapse = "|")

group_mapping <- list(
  `Pesticidal crystal` = "pest_crys",
  `Neurotoxin` = c("scoloptoxin", "latrotoxin", "conotoxin"),
  `Metalloprotease` = c("metalloproteinase", "metallopeptidase", "metalloendopeptidase", "serralysin"),
  `Porin` = "neoverrucotoxin"
)
TOXIN_GROUP_MAP <- lapply(names(group_mapping), \(x) lapply(group_mapping[[x]], \(y) setNames(x, y))) |> unlist()

resolve_toxins <- function(row) {
  h <- row$header
  if ((str_detect(h, "toxin") && str_detect(h, "Chironex")) ||
    row$toxin_group == "cftx") {
    return("Jft")
  }
  if (row$toxin_group %in% names(TOXIN_GROUP_MAP)) {
    return(TOXIN_GROUP_MAP[row$toxin_group])
  }
  if (row$toxin_group != "toxin") {
    return(str_to_title(row$toxin_group))
  }
  "Unspecified toxin"
}


# TODO: Redo this entirely match PFAM and PANTHER using accessions
# but use keyword searches for header and eggNOG description
toxin_results <- M$data |>
  filter(
    if_any(any_of(ANNO_COLS), \(x) grepl(toxin_query, x, ignore.case = TRUE))
    # if_any(contains("GO"), \(x) x == "toxin_related")
  ) %>%
  group_by(GroupUP) |>
  summarise(across(everything(), \(x) paste0(x, collapse = ";"))) %>%
  bind_cols(
    select(., all_of(ANNO_COLS)) |> unite(descriptors, sep = ";")
  ) |>
  select(
    ProteinId, header, all_of(ANNO_COLS), contains("GO"), pcoverage_align, descriptors
  ) |>
  mutate(
    descriptors = lapply(
      descriptors,
      \(x) str_split_1(x, ";") |> discard(\(x) x == "NA")
    ),
    toxin_group = map_chr(descriptors, \(x) grep_most_frequent(toxin_words, x)),
  ) |>
  filter(!grepl("botulinum toxin substrate", header) &
    !grepl("botulinum toxin substrate", eggNOG_description))


toxin_results$toxin_group <- apply(toxin_results, 1, resolve_toxins)


toxin_groups <- table(toxin_results$toxin_group) |>
  table2tb(id_col = "Group") |>
  mutate(percentage = n / sum(n))

toxin_pie <- plot_ly(toxin_groups,
  labels = ~Group, values = ~percentage, type = "pie",
  textinfo = "label+percent"
) |> layout(
  title = "Toxin group proportions", xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
)


toxin_waffle <- ggplot(toxin_groups, aes(values = n, fill = Group)) +
  scale_fill_paletteer_d(PALETTE) +
  geom_waffle(n_rows = 10, size = 0.5, make_proportional = FALSE) +
  guides(fill = guide_legend("Toxin group"))


save(c(TABLES, GRAPHS), glue("{M$outdir}/Figures/taxonomy"))
