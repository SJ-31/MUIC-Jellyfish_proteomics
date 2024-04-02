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


run <- runData("./results/C_indra", TRUE)
GRAPHS <- list()
TABLES <- list()


# Coverage metrics
cov_align <- compareFirstSecL(run, "coverage_alignlen",
                              TRUE, "ProteinId")
num_peps <- compareFirstSecW(run, "num_peps", "ProteinId", TRUE)
test_num_peps <- wilcoxWrapper(num_peps)
test_coverage <- wilcoxWrapper(compareFirstSecW(run, "coverage_alignlen", "ProteinId"))
GRAPHS$run_coverage <- passDensityPlot(cov_align, 0.05) + labs(x = "percent coverage")


# Unique proteins to each run
run_uniques <- passUniques(run)
percent_found <- dplyr::bind_cols(notMissing(run$first),
                                  notMissing(run$sec)) %>%
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
# Test if modifications are associated with category???

# Annotation metrics
counts <- list()
counts$first <- getCounts(run$first)
counts$sec <- getCounts(run$sec)
# sv_first <- counts$first$go$GO %>%
#   purrr::map_dbl(getSV)
# sv_sec <- counts$sec$go$GO %>% purrr::map_dbl(getSV)
# sv_compare <- tibble(sv = sv_first, pass = "first") %>%
#   bind_rows(tibble(sv = sv_sec, pass = "sec"))

wanted_cols <- c("GO_counts", "GO_max_sv", "num_peps")

avgStdevs <- function(tb, cols, stat) {
  result <- vector()
  for (col in cols) {
    result <- c(result, stat(tb[[col]]))
  }
  return(result)
}

avgStdevs(run$first, wanted_cols, base::mean)

per_protein <- tibble(
  metric = c("Stdev GO count", "Average GO count", "Stdev semantic value", "Avg (max) semantic value", "Stdev peptide count", "Avg peptide count")
)

GO_counts <- list()


# Counts of go terms per protein were recorded already in combine_all.r

# Max sv per protein (better do this in combine_all

p_rename <- c(NCBI_ID = "Accession Number", coverage_nmatch.prev = "Sequence coverage [%]")

# Compare with previous
p_all <- read_tsv("./data/reference/previous_all.tsv") %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))
p_toxins <- read_tsv("./data/reference/previous_toxins.tsv") %>%
  rename(., all_of(p_rename)) %>%
  select(-contains(" "))

compare_all <- inner_join(p_all, run$first, by = join_by(NCBI_ID)) %>%
  left_join(., run$sec, by = join_by(NCBI_ID), suffix = JOIN_SUFFIX) %>%
  select(c(NCBI_ID, header.first, category.first, coverage_nmatch.prev, coverage_nmatch.first, coverage_nmatch.sec)) %>%
  mutate(coverage_nmatch.first = coverage_nmatch.first * 100,
         coverage_nmatch.sec = coverage_nmatch.sec * 100)

cov_longer <- compare_all %>%
  select(contains("coverage"), NCBI_ID) %>%
  rename_with(., \(x) purrr::map_chr(x, \(y) {
    if (grepl("NCBI_ID", y)) return(y)
    y <- gsub("coverage_nmatch", "", y)
    if (y == ".prev") return("previous")
    else if (y == ".first") return("first")
    return("second")
  })) %>%
  pivot_longer(cols = !NCBI_ID) %>%
  mutate(value = round(value, 2))

GRAPHS$shared_cov_bp <- cov_longer %>% ggplot(aes(y = value, fill = name)) +
  geom_boxplot() +
  theme(axis.text.x = element_blank()) +
  xlab("Seqence coverage") +
  ggtitle("Shared protein sequence coverage")

# For space purposes, drop entries where ALL had less than 30% sequence coverage, as well as any entries both second and first passes cannot be compared
threshold <- 0.4
cov_longer <- cov_longer %>%
  group_by(NCBI_ID) %>%
  nest()
mask <- cov_longer %>% apply(1, \(x) {
  data <- x$data
  if (all(data$value < threshold) | any(is.na(data$value))) return(FALSE)
  return(TRUE)
})
cov_longer <- cov_longer[mask,] %>% unnest()
cov_longer$name <- factor(cov_longer$name, levels = c("previous", "first", "second"))

GRAPHS$protein_wise_coverage <- cov_longer %>%
  ggplot(aes(x = NCBI_ID, y = value, fill = name)) +
  geom_bar(stat = "identity") +
  stat_identity(geom = "text", colour = "black", size = 5, aes(label = value), position = position_stack(vjust = 0.5)) +
  ylab("Sequence coverage") +
  xlab(element_blank()) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# Filter out new proteins
compare_toxin <- compare_all %>% filter(NCBI_ID %in% p_toxins$NCBI_ID)
new_proteins <- run$first %>% filter(!NCBI_ID %in% compare_all$NCBI_ID)
new_toxins <- new_proteins %>% filter(category == "venom_component")



