if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
SECOND_PASS_ENGINES <- c("identipy", "msgf", "msfragger", "comet")

TABLES <- list()

db <- glue("{M$chosen_path}/Databases/decoysWnormal.fasta")
db_size <- function(fasta_file) {
  seqkit_stat(fasta_file)$num_seqs[[1]] %>% as.numeric()
}

complete_db_size <- db_size(db)
sizes <- lapply(SECOND_PASS_ENGINES, \(x) {
  db_size(glue("{M$chosen_path}/2-Second_pass/BK_databases/{x}_bk_database.fasta"))
}) %>% `names<-`(SECOND_PASS_ENGINES)

get_percolator <- function(pass) {
  dir <- glue("{M$chosen_path}/{pass}/Percolator")
  lapply(SECOND_PASS_ENGINES, \(x) {
    read_tsv(glue("{dir}/{x}_percolator_proteins.tsv"))
  }) %>% `names<-`(SECOND_PASS_ENGINES)
}

first_pass <- get_percolator("1-First_pass")
sec_pass <- get_percolator("2-Second_pass")

get_stats <- function(engine_list, pass) {
  get_row <- function(lst_subset) {
    tb <- lst_subset[[1]]
    tibble(
      engine = names(lst_subset),
      num_proteins = nrow(tb),
      num_passed_fdr = tb %>% filter(`q-value` < M$fdr) %>% nrow(),
      num_passed_pep = tb %>% filter(`posterior_error_prob` < M$fdr) %>% nrow()
    )
  }
  lmap(engine_list, get_row) %>%
    bind_rows() %>%
    mutate(pass = pass)
}


stats <- bind_rows(get_stats(first_pass, "first"), get_stats(sec_pass, "sec"))
SEQ_MAP <- read_tsv(glue("{M$chosen_path}/Databases/seq-header_mappings.tsv"))

NEW_IN_SEC <- c()

ALL_RESULTS <- bind_rows(
  mutate(M$run$first, pass = "first"),
  mutate(M$run$second, pass = "sec")
) %>%
  separate_longer_delim(., "MatchedPeptideIds", ";") |>
  inner_join(SEQ_MAP,
    by = join_by(x$MatchedPeptideIds == y$id),
    suffix = c("", ".y")
  ) |>
  rename(MatchedPeptideIdsHeader = header.y) |>
  select(-contains(".y"))


KEPT_PROTEINS <- c(ALL_RESULTS$header, ALL_RESULTS$MatchedPeptideIdsHeader) |> unique()
ALL_RESULTS <- local({
  expanded <- ALL_RESULTS %>%
    mutate(header = MatchedPeptideIdsHeader) %>%
    filter(!is.na(header))
  bind_rows(ALL_RESULTS, expanded)
}) %>%
  distinct(header, .keep_all = TRUE)

# Per-engine analysis function
per_engine <- function(engine) {
  engine_results <- list()
  data <- dplyr::bind_rows(
    separate_longer_delim(first_pass[[engine]], "ProteinId", ",") %>%
      mutate(., pass = "first"),
    separate_longer_delim(sec_pass[[engine]], "ProteinId", ",") %>%
      mutate(., pass = "sec"),
  ) %>%
    inner_join(., SEQ_MAP, by = join_by(x$ProteinId == y$id))

  data <- data %>% mutate(
    from = dplyr::case_when(
      str_detect(ProteinId, "P") ~ "database",
      str_detect(ProteinId, "D") ~ "denovo",
      str_detect(ProteinId, "T") ~ "transcriptome",
    ),
    num_peps = map_dbl(peptideIds, \(x) str_count(x, " ") + 1)
  )

  first <- data %>% filter(pass == "first")
  sec <- data %>% filter(pass == "sec")
  data <- mutate(data, lost_in_sec = !header %in% sec$header)


  lost_in_sec <- filter(first, !ProteinId %in% sec$ProteinId)
  # print(all(sec$header %in% first$header))
  NEW_IN_SEC <<- c(NEW_IN_SEC, filter(sec, !header %in% first$header) |> pluck("ProteinId")) # Should be empty

  found_in_both <- inner_join(first, sec,
    by = join_by(x$header == y$header),
    suffix = c(".first", ".sec")
  ) %>%
    select(-contains("pass")) %>%
    select(-from.first) %>%
    rename(from = from.sec)

  found_in_both <- found_in_both %>% mutate(
    delta_num_peps = num_peps.sec - num_peps.first,
    delta_q_value = `q-value.sec` - `q-value.first`,
    delta_pep = posterior_error_prob.sec - posterior_error_prob.first
  )
  # Want to check the number of peptides in the proteins lost

  # Run pairwise tests and get summary statistics to check
  # if the second pass performed better than the first
  to_test <- c("num_peps", "q-value", "posterior_error_prob")
  engine_results$pairwise_tests <- pairwise_tests_tb(
    found_in_both, to_test, c("less", "greater", "greater"),
    \(x, y, ...) wilcox.test(x, y, paired = TRUE, ...)
  )

  tab <- table(data$from, data$pass) %>% table2df()
  chi <- chisq.test(tab)
  # Check if the distribution of protein types identified differs between
  # runs. It shouldn't, because the second should only identify proteins in the
  # first
  chi$data.name <- "Frequency of proteins from different sources"

  chi <- htest2tb(chi)

  kept_data <- data %>% filter(header %in% KEPT_PROTEINS)
  merged <- inner_join(kept_data, ALL_RESULTS, by = join_by(header)) %>%
    select(-contains(".x|.y"))

  sources <- c("transcriptome", "denovo", "database")
  engine_results$OR <- lapply(
    sources,
    \(x) {
      data$lost_in_sec
      table <- table(data$from != x, data$lost_in_sec)
      # print(table)
      or <- table %>% get_odds_ratio()
      upper <- table %>% get_odds_ratio(CI = TRUE)
      lower <- table %>% get_odds_ratio(CI = TRUE, side = "lower")
      tibble(source = x, odds_ratio = or, CI_lower = lower, CI_upper = upper)
    }
  ) %>%
    bind_rows()
  # Setup means that we interpret OR as the odds of not being lost in the
  # second pass is OR times as high in proteins of the current group than
  # all others

  # Change to check peptides
  # Check if the lengths of proteins lost have a statistically significant difference
  # with those kept
  as_peps <- data |>
    mutate(peptideIds = map_chr(peptideIds, fill_peptide_gaps)) |>
    separate_longer_delim("peptideIds", ";") |>
    mutate(length = nchar(peptideIds))
  lengths_lost <- data %>%
    filter(lost_in_sec) %>%
    purrr::pluck("length")
  lengths_kept <- data %>%
    filter(!lost_in_sec) %>%
    purrr::pluck("length")

  test <- wilcox.test(lengths_lost, lengths_kept)
  test$data.name <- "lost x retained"
  p <- test$`p.value`
  test <- htest2tb(test)
  if (p < 0.05) {
    t2 <- wilcox.test(lengths_lost, lengths_kept, alternative = "less")
    t2$data.name <- "lost x retained"
    t2$alternative <- "lost less"
    test <- bind_rows(test, htest2tb(t2))
  }
  engine_results$htests <- bind_rows(chi, test)
  return(engine_results)
}


TABLES <- purrr::reduce(SECOND_PASS_ENGINES, \(acc, x) {
  mutate_bind <- function(col) {
    acc[[col]] <- bind_rows(
      acc[[col]],
      mutate(results[[col]], engine = x)
    )
  }
  results <- per_engine(x)
  acc$pairwise_tests <- mutate_bind("pairwise_tests")
  acc$htests <- mutate_bind("htests")
  acc$OR <- mutate_bind("OR")
  return(acc)
}, .init = list(pairwise_tests = tibble(), htests = tibble(), OR = tibble()))

pw <- TABLES$pairwise_tests
pw$conclusion <- pmap(
  list(pw$alternative, pw$two_sided_significant, pw$alternative_significant),
  \(alt, two_sided_significant, alternative_significant) {
    alt <- ifelse(str_detect(alt, "less"), "less", "greater")
    if (two_sided_significant != "Y" || is.na(two_sided_significant)) {
      NA
    } else if (alternative_significant == "Y") {
      glue("first {alt}")
    } else {
      glue("second {alt}")
    }
  }
) |> unlist()

TABLES$pairwise_tests_raw <- pw

sub <- substitute_all(
  c("num_peps", "pcoverage_align"),
  c("peptide number", "percent coverage"),
  \(x) gsub("_", " ", x)
)

TABLES$pairwise_tests <- pw %>%
  filter(two_sided_significant == "Y") |>
  filter(metric == "num_peps") |>
  mutate(alternative = map_chr(alternative, \(x) str_remove(x, ".first"))) |>
  rename(significant = alternative_significant) |>
  select(-c(mean_diff, two_sided_p_value, two_sided_significant)) |>
  rename_with(
    \(vec) {
      rename_helper <- function(x) {
        if (str_detect(x, "sided")) {
          x <- str_replace(x, ".sided", "-sided")
        }
        if (str_detect(x, "p_value")) {
          x <- str_replace(x, "p_value", "p-value")
        }
        str_replace_all(x, "_", " ")
      }
      map_chr(vec, rename_helper)
    },
    everything()
  ) |>
  gt(rowname_col = "metric") %>%
  fmt(., columns = "metric", fns = \(x) map_chr(x, sub)) %>%
  text_case_match(
    "two.sided" ~ "two-sided",
    ".first" ~ "first",
    "_" ~ " pass ",
    .replace = "partial"
  ) %>%
  fmt_number(., columns = contains("p-value"), decimals = 5) %>%
  tab_stubhead(label = "Metric")


TABLES$stats <- stats

# ----------------------------------------
# Confirm that new proteins were seen previously
get_blast <- function(pass, prefix = M$chosen_prefix, path = M$chosen_path) {
  bind_rows(
    read_tsv(glue("{path}/{pass}/Unmatched/BLAST/{prefix}_blast_matched.tsv")),
    read_tsv(glue("{path}/{pass}/Unmatched/BLAST/{prefix}_blast_unmatched.tsv"))
  )
}

new_in_sec <- M$run$second |>
  filter(!header %in% M$run$first$header) |>
  distinct(ProteinId, .keep_all = TRUE) |>
  pluck("ProteinId")
bmatches <- get_blast("1-First_pass")
bmatches2 <- get_blast("2-Second_pass")

b_ids <- c(bmatches$ProteinId, flatten_by(bmatches$MatchedPeptideIds, ";"))
b_ids2 <- c(bmatches2$ProteinId, flatten_by(bmatches2$MatchedPeptideIds, ";"))
found_previously <- new_in_sec %>% keep(\(x) x %in% b_ids)
found_in_blast <- new_in_sec %>% keep(\(x) (x %in% b_ids2) && !(x %in% found_previously))
others <- new_in_sec |> discard(\(x) x %in% found_previously || x %in% found_in_blast)

TABLES$new_in_second_pass <- glue("
Total number of new proteins in second pass: {length(b_ids2)}
Number found in first pass: {length(found_previously)}
% found in first pass: {length(found_previously)/length(b_ids)}
Number found from blast in second pass {length(found_in_blast)}
% found from blast in second pass {length(found_in_blast)/length(b_ids2)}
")

TABLES$peptide_length_tests <- TABLES$htests |>
  get_adjusted_p() |>
  filter(!grepl("Chi", method)) |>
  rename(pair = data) |>
  mutate(pair = paste0(pair, " (", engine, ")")) |>
  conclude_one_sided() |>
  pairwise_conclusion2gt()

GRAPHS <- list()
GRAPHS$odds_ratio <- TABLES$OR |> ggplot(aes(x = engine, y = odds_ratio, fill = source)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, position = position_dodge(.9)) +
  scale_fill_paletteer_d("awtools::ppalette") +
  ylab("Odds of being lost in second pass") +
  guides(fill = guide_legend("Protein type"))
attr(GRAPHS$odds_ratio, "height") <- 6

save(c(GRAPHS, TABLES), glue("{M$outdir}/pass_differences"))
