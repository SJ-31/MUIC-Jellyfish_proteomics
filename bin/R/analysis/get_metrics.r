library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)
args <- list(r_source = "./bin/R")
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/GO_helpers.r"))
source(glue("{args$r_source}/analysis/metric_functions.r"))


EGGNOG_COLS <- c("EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")


run <- runData("C_indra", TRUE, "./results/C_indra_A")
GRAPHS <- list()
TABLES <- list()


# Coverage metrics
cov_align <- compareFirstSecL(
  run, "pcoverage_align",
  TRUE, "ProteinId"
)
num_peps <- compareFirstSecW(run, "num_peps", "ProteinId", TRUE)
GRAPHS$run_coverage <- passDensityPlot(cov_align, 0.05) + labs(x = "percent coverage")


# Unique proteins to each run
run_uniques <- passUniques(run)
percent_found <- dplyr::bind_cols(
  notMissing(run$first),
  notMissing(run$sec)
) %>%
  `colnames<-`(c("first", "sec")) %>%
  tibble::rownames_to_column(., var = "metric") %>%
  as_tibble()
wanted <- c("lineage", "Mods", "flashlfq_mean", "maxlfq_mean", "CAZy", "PFAMs", "EC", "BRITE", "interpro_accession", "PANTHER", "eggNOG_OGs", "UniProtKB_ID", "GO", "directlfq_mean")
GRAPHS$percent_found <- percent_found %>%
  dplyr::filter(!(first == 100 & sec == 100) & metric %in% wanted) %>%
  pivot_longer(cols = c("first", "sec")) %>%
  ggplot(aes(x = metric, y = value, fill = name)) +
  geom_bar(position = "dodge", stat = "identity") +
  ylab("% not missing")

# --------------------------------------------------------

# PTMs
has_mods <- run$first %>% filter(!is.na(Mods))
UNIQUE_MODS <- run$first$Mods %>%
  discard(is.na) %>%
  map(\(x) str_split_1(x, "\\|")) %>%
  unlist() %>%
  map_chr(\(x) gsub(" [1-9]+", "", x)) %>%
  unique() %>%
  discard(\(x) grepl("0$", x))

ptms_first <- modMetrics(run$first)
ptms_sec <- modMetrics(run$sec)
ptm_percent_diff <- abs(ptms_sec$percentages - ptms_first$percentages)

# Test if modifications are associated using chi square

# Transform tb so that proteins are in rows and mods in cols
# A protein has TRUE in the mod col if it has that mod
tb <- run$first
chi <- tb %>%
  filter(category != "other") %>%
  select(ProteinId, category) %>%
  left_join(., df2Tb(ptms_first$count_df, "ProteinId")) %>%
  mutate(across(is.double, \(x) ifelse(is.na(x), FALSE, TRUE)))

no_mods <- chi %>%
  select(-c(ProteinId, category)) %>%
  apply(1, \(x) ifelse(any(x), FALSE, TRUE)) %>%
  unlist()
chi$none <- no_mods

# Create 2 x 2 contingency table with mods as columns and categories as rows
# But since mods are not mutually exclusive, will need to test each
# modification individually against the categories
chosen <- "Met_Oxidation"
tab <- table(chi$category, chi[[chosen]])
ptm_tests <- list()
for (mod in UNIQUE_MODS) {
  print(glue("Testing {mod}"))
  ptm_tests[[mod]]$test <- chisq.test(chi$category, chi[[mod]])
  ptm_tests[[mod]]$table <- table(chi$category, chi[[mod]])
  print(ptm_tests[[mod]])
}


# --------------------------------------------------------


# Annotation metrics
counts <- list()
counts$first <- getCounts(run$first)
counts$sec <- getCounts(run$sec)
# sv_first <- counts$first$go$GO %>%
#   purrr::map_dbl(getSV)
# sv_sec <- counts$sec$go$GO %>% purrr::map_dbl(getSV)
# sv_compare <- tibble(sv = sv_first, pass = "first") %>%
#   bind_rows(tibble(sv = sv_sec, pass = "sec"))

wanted_cols <- c("GO_counts", "GO_max_sv", "num_peps", "pcoverage_align")

avgStdevs <- function(tb, cols, stat) {
  result <- vector()
  for (col in cols) {
    result <- c(result, stat(tb[[col]]))
  }
  return(result)
}


test_num_peps <- wilcoxWrapper(num_peps)
test_coverage <- wilcoxWrapper(compareFirstSecW(run, "pcoverage_align", "ProteinId"))
# Results per protein
# Evaluate significance of each
per_protein <- tibble(
  metric = rep(wanted_cols, 2),
  type = c(rep("mean", length(wanted_cols)), rep("stdev", length(wanted_cols))),
  first = c(
    avgStdevs(run$first, wanted_cols, \(x) mean(x, na.rm = TRUE)),
    avgStdevs(run$first, wanted_cols, \(x) sd(x, na.rm = TRUE))
  ),
  sec = c(
    avgStdevs(run$sec, wanted_cols, \(x) mean(x, na.rm = TRUE)),
    avgStdevs(run$sec, wanted_cols, \(x) sd(x, na.rm = TRUE))
  ),
  wilcox_result = c() # TODO: Add in this column after the results have
  # come in
)
# per_protein %>%
#   mutate(metric = paste0(metric, "_", type)) %>%
#   rename()
# pivot_longer(cols = c(first, sec)) %>%
#   ggplot(aes(x = metric, y = value, fill = name)) +
#   geom_bar(stat = "identity", position = "dodge")


# Compare with previous results
p_rename <- c(NCBI_ID = "Accession Number", pcoverage_nmatch.prev = "Sequence coverage [%]")
p_all <- read_tsv("./data/reference/previous_all.tsv") %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))
p_toxins <- read_tsv("./data/reference/previous_toxins.tsv") %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))

compare_all <- inner_join(p_all, run$first, by = join_by(NCBI_ID)) %>%
  left_join(., run$sec, by = join_by(NCBI_ID), suffix = JOIN_SUFFIX) %>%
  select(c(NCBI_ID, header.first, category.first, pcoverage_nmatch.prev, pcoverage_nmatch.first, pcoverage_nmatch.sec)) %>%
  mutate(
    pcoverage_nmatch.first = pcoverage_nmatch.first * 100,
    pcoverage_nmatch.sec = pcoverage_nmatch.sec * 100
  )

cov_longer <- compare_all %>%
  select(contains("coverage"), NCBI_ID) %>%
  rename_with(., \(x) purrr::map_chr(x, \(y) {
    if (grepl("NCBI_ID", y)) {
      return(y)
    }
    y <- gsub("pcoverage_nmatch", "", y)
    if (y == ".prev") {
      return("previous")
    } else if (y == ".first") {
      return("first")
    }
    return("second")
  })) %>%
  mutate(first_diff = first - previous,
         second_diff = second - previous) %>%
  select(-c(first, second, previous)) %>%
  pivot_longer(cols = !NCBI_ID) %>%
  mutate(value = round(value, 2))

GRAPHS$shared_cov_bp <- cov_longer %>% ggplot(aes(y = value, fill = name)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  xlab("Seqence coverage") +
  ggtitle("Shared protein sequence coverage")

# For space purposes, drop entries where ALL had less than 30% sequence coverage, as well as any entries both second and first passes cannot be compared
threshold <- 0
cov_longer <- cov_longer %>%
  group_by(NCBI_ID) %>%
  nest()
mask <- cov_longer %>% apply(1, \(x) {
  data <- x$data
  if (all(data$value < threshold) | any(is.na(data$value))) {
    return(FALSE)
  }
  return(TRUE)
})
cov_longer <- cov_longer[mask,] %>% unnest()
cov_longer$name <- factor(cov_longer$name, levels = c("previous", "first", "second"))

GRAPHS$protein_wise_coverage <- cov_longer %>%
  ggplot(aes(x = NCBI_ID, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Sequence coverage change from previous (%)") +
  xlab(element_blank()) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )


# Filter out new proteins
compare_toxin <- compare_all %>% filter(NCBI_ID %in% p_toxins$NCBI_ID)
new_proteins <- run$first %>% filter(!NCBI_ID %in% compare_all$NCBI_ID)
new_toxins <- new_proteins %>% filter(category == "venom_component")


# --------------------------------------------------------
#' Split and a ProteinGroupId string by the ";", optionally remove the
#' numbers and leave unique groups
#'
splitGroupStr <- function(group_str, remove_nums = FALSE, unique = FALSE) {
  if (remove_nums)
    group_str <- gsub("[0-9]+", "", group_str);
  split <- str_split_1(group_str, ";")
  if (unique) split <- base::unique(split)
  return(split)
}

# Engine analysis
num_ids <- bind_rows(run$first, run$sec) %>%
  filter(ProteinGroupId != "U") %>%
  select(c(ProteinGroupId, pcoverage_nmatch, ProteinId, num_peps, num_unique_peps)) %>%
  mutate(engine_count = purrr::map_dbl(ProteinGroupId, \(x) {
    x <- splitGroupStr(x, TRUE, TRUE) %>%
      discard(\(x) x == "U")
    return(length(x))
  }))

# Figure out which engines had the biggest contributions
engine_counts <- num_ids$ProteinGroupId %>%
  lapply(., \(x) splitGroupStr(x, TRUE)) %>%
  unlist() %>%
  discard(\(x) x == "U") %>%
  table()

# Is there an association between the number of matched peptides and the
# identity of the engines matching them? That is, are some engines more
# "isolated" than others? i.e. engine B tends to identify proteins that other
# engine don't. But also remember that in your setup each protein needs to
# be identified by at least two standard engines (only one in open search)
tb <- run$first
engine_hits <- tb %>%
  select(ProteinId, ProteinGroupId) %>%
  apply(1, \(x) {
    groups <- splitGroupStr(x["ProteinGroupId"],
                            TRUE, TRUE)
    row <- tibble(ProteinId = x["ProteinId"])
    for (g in groups) {
      row <- add_column(row, !!g := TRUE)
    }
    return(row)
  }) %>%
  bind_rows() %>%
  replaceNaAll(., FALSE)
# Collapse the contingency table?
engine_hits <- engine_hits %>%
  inner_join(., dplyr::select(tb, category, ProteinId)) %>%
  select(-ProteinId) %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  filter(n() > 30) %>%
  ungroup()


categories <- unique(engine_hits$category)
current <- "comet"
tl <- table(engine_hits$category, engine_hits[[current]])
test <- chisq.test(tl)
# If there is some significant association, collapse the table
# in order to find which one specifically
if (test$p.value < 0.05) {
  for (cat in categories) {
    t <- table(engine_hits$category == cat, engine_hits[[current]])
    tst <- chisq.test(t)
    print(tst)
  }
}
# TODO: Wrap this in a function


# Correlation between no. identifications by engines and coverage
engine_cor <- cor.test(num_ids$engine_count, num_ids$pcoverage_nmatch) # A weak positive correlation, but statistically significant
# Correlation between peptide number and coverage
n_peps_cor <- cor.test(num_ids$num_unique_peps, num_ids$pcoverage_nmatch) # A weak positive correlation, but statistically significant


# --------------------------------------------------------

minMaxScaler <- function(vec) {
  max <- max(vec, na.rm = TRUE)
  min <- min(vec, na.rm = TRUE)
  return(purrr::map_dbl(vec, \(x) (x - min) / (max - min)))
}

