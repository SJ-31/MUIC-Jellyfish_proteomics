allIntervals <- function(ids) {
  getIdIntervals <- function(id) {
    lapply(ENGINES, \(x) py$tracer$engineAlignmentIntervals(id, x)) %>%
      `names<-`(ENGINES)
  }
  lapply(ids, getIdIntervals) %>% `names<-`(ids)
}

intervalCoverage <- function(interval_list, length, id) {
  sumInterval <- function(intervals) {
    if (length(intervals) == 0) {
      return(0)
    }
    s <- sum(map_dbl(intervals, \(x) x[2] - x[1]))
    s / length
  }
  lmap(
    interval_list,
    \(x) tibble(engine = names(x), ProteinId = id, coverage = sumInterval(x[[1]]))
  ) %>%
    bind_rows() %>%
    pivot_wider(names_from = engine, values_from = coverage)
}

# File and path setup
pass <- "1-First_pass"
PERCOLATOR_DIR <- c(glue("{PATH}/{pass}/Percolator"), glue("{PATH}/{pass}/Open_search/Percolator"))
open_search_engines <- c("metamorpheusGTPMD", "msfraggerGPTMD", "msfraggerGlyco")
percolator_files <- list.files(PERCOLATOR_DIR,
  pattern = "_percolator_proteins.tsv", full.names = TRUE
)
alignments <- alignmentData(PATH, "first")
combined_results <- run$first

# For TESTING
combined_results <- combined_results %>% slice(1:3000)


# Engine distribution analysis
ENGINES <- percolator_files %>%
  map_chr(., \(x) {
    gsub("_percolator_proteins.tsv", "", x) %>%
      gsub(".*/", "", .)
  })

names(percolator_files) <- ENGINES
percolator_files <- as.list(percolator_files)
denovo_transcriptome_ids <- combined_results$MatchedPeptideIds %>%
  discard(is.na) %>%
  lapply(., str_split_1, pattern = ";") %>%
  unlist()
database_ids <- combined_results$ProteinId
all_found_protein_ids <- c(database_ids, denovo_transcriptome_ids)
standard_search <- combined_results %>% filter(ID_method == "standard")

engine_stats <- tibble()
percolator_tibbles <- local({
  getEngine <- function(file) {
    tb <- tibbleDuplicateAt(read_tsv(file[[1]]), "ProteinId", ",") %>% mutate(engine = names(file))
    filtered <- dplyr::filter(tb, ProteinId %in% all_found_protein_ids) %>%
      mutate(num_peps = map_dbl(peptideIds, \(x) length(str_split_1(x, " "))))
    row <- tibble(
      engine = names(file),
      n_groups_all = length(unique(tb$ProteinGroupId)),
      n_groups_accepted = length(unique(filtered$ProteinGroupId)),
      n_proteins_all = nrow(tb),
      n_proteins_below_fdr = tb %>% filter(`q-value` <= FDR) %>% nrow(),
      n_proteins_accepted = nrow(filtered),
      # To be accepted, a protein needs to be matched by at least two other engines
      # And the q value of at least two psm matches needs to be below the FDR
      n_other_proteins = map_lgl(tb$ProteinId, \(x) str_detect(x, "D|T")) %>% sum(),
    ) %>%
      mutate(n_database = n_proteins_all - n_other_proteins)
    engine_stats <<- dplyr::bind_rows(engine_stats, row)
    return(filtered)
  }
  lmap(percolator_files, getEngine) %>% `names<-`(ENGINES)
})

normal_engine_tbs <- percolator_tibbles[!names(percolator_tibbles) %in% open_search_engines]
percolator_tb <- dplyr::bind_rows(percolator_tibbles)

num_peptides_matched <- apply(combined_results, 1, \(x) {
  protein_id <- x["ProteinId"]
  ids_to_match <- protein_id
  matched_peptides <- x["MatchedPeptideIds"]
  if (!is.na(matched_peptides)) {
    ids_to_match <- c(ids_to_match, str_split_1(matched_peptides, ";"))
  }
  row <- tibble(ProteinId = protein_id)
  for (e in ENGINES) {
    current <- percolator_tibbles[[e]] %>% filter(ProteinId %in% ids_to_match)
    if (nrow(current) != 0) {
      col <- tibble(!!e := sum(current$num_peps))
    } else {
      col <- tibble(!!e := 0)
    }
    row <- bind_cols(row, col)
  }
  return(row)
}) %>%
  bind_rows()

# Determine which engines contributed the most to protein coverage
# Do so by multivariate linear regression -> it makes sense that
with_normal <- local({
  normal <- num_peptides_matched %>% select(all_of(c("ProteinId", names(normal_engine_tbs))))
  standard_search %>%
    select(., c(ProteinId, pcoverage_align)) %>%
    inner_join(., normal, by = join_by("ProteinId"))
})

reticulate::source_python(glue("{args$python_source}/trace_alignments.py"))
matched_peptides <- select(combined_results, c(ProteinId, MatchedPeptideIds))
protein_tb <- normal_engine_tbs %>%
  bind_rows() %>%
  tibbleDuplicateAt(., "peptideIds", " ") %>%
  mutate(peptideIds = map_chr(peptideIds, cleanPeptide))

py$tracer <- AlignmentTracer(protein_tb, alignments$peptides, matched_peptides)
# Try to measure the proportion of alignments each engine takes up with its peptides
traced_intervals <- allIntervals(combined_results$ProteinId)
interval_coverage <- pmap(
  list(traced_intervals, combined_results$length, combined_results$ProteinId),
  \(x, y, z) intervalCoverage(x, y, z)
) %>% bind_rows()
interval_groups <- pmap(
  list(traced_intervals, combined_results$ProteinId),
  \(x, y) getIntervalGroups(x, y)
) %>%
  bind_rows() %>%
  as_tibble()

standard_engine_coverage <- interval_coverage %>%
  filter(ProteinId %in% standard_search$ProteinId) %>%
  select(-all_of(open_search_engines))

GRAPHS$engine_peptide_coverage <- standard_engine_coverage %>%
  pivot_longer(., cols = -ProteinId) %>%
  ggplot(aes(y = value, x = name, fill = name)) +
  geom_boxplot() +
  labs(title = "Density of per-engine peptide coverage") +
  xlab("Coverage (%)") +
  ylab("Engine name") +
  theme(legend.position = "none")

# Identify which engine, if any, is the best-performing
cov_list <- standard_engine_coverage %>%
  select(-ProteinId) %>%
  as.list()
ks <- kruskal.test(cov_list)
combos <- combn(names(cov_list), 2)
tests <- lapply(
  seq_len(ncol(combos)),
  \(x) wilcox.test(cov_list[[combos[1, x]]], cov_list[[combos[2, x]]])
)
test_tb <- tibble(
  p_value = map_dbl(tests, \(x) x$p.value),
  p_adjust = p.adjust(p_value),
  pair = map_chr(
    seq_len(ncol(combos)),
    \(x) paste0(combos[1, x], " x ", combos[2, x])
  ),
  significant = map_dbl(p_adjust, \(x) x < 0.05)
)

# id_list <- lapply(engines, function(x)
#   tibbles %>% bind_rows() %>% filter(engine == x))$ProteinId
# }) %>% `names<-`(engines)
# svg("engine_intersections.svg")
# venn(id_list,
#   plotsize = 100, borders = FALSE, box = FALSE,
#   ilcs = 1, sncs = 1.5,
#   zcolor = "red, blue, green, white, orange, yellow"
# )

# filtered_tibbles <- lmap(tibbles, \(x) inner_join(x[[1]], final_results, by = join_by("ProteinId")))

# TODO: Extract different columns (e.g. lineage, GO) from "filtered tibbles" to see
# if there are statistically significant differences between the different engines
# Are the engines biased in how they find proteins??
