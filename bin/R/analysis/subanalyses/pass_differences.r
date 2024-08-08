# BUG Don't really need this one tbh
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
SECOND_PASS_ENGINES <- c("identipy", "msgf", "msfragger", "comet")

TABLES <- list()

db <- glue("{M$path}/Databases/decoysWnormal.fasta")
db_size <- function(fasta_file) {
  seqkit_stat(fasta_file)$num_seqs[[1]] %>% as.numeric()
}

complete_db_size <- db_size(db)
sizes <- lapply(SECOND_PASS_ENGINES, \(x) {
  db_size(glue("{M$path}/2-Second_pass/BK_databases/{x}_bk_database.fasta"))
}) %>% `names<-`(SECOND_PASS_ENGINES)

get_percolator <- function(pass) {
  dir <- glue("{M$path}/{pass}/Percolator")
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
SEQ_MAP <- read_tsv(glue("{M$path}/Databases/seq-header_mappings.tsv"))


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


  new_in_sec <- sec %>% filter(!header %in% first$header) # Should be empty

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

  # Check if the lengths of proteins lost have a statistically significant difference
  # with those kept
  lengths_lost <- data %>%
    filter(lost_in_sec) %>%
    purrr::pluck("length")
  lengths_kept <- data %>%
    filter(!lost_in_sec) %>%
    purrr::pluck("length")

  test <- wilcox.test(lengths_lost, lengths_kept)
  test$data.name <- "Lengths of proteins lost in second pass vs lengths of those retained"
  p <- test$`p.value`
  test <- htest2tb(test)
  if (p < 0.05) {
    t2 <- wilcox.test(lengths_lost, lengths_kept, alternative = "less")
    t2$data.name <- "Lengths of proteins lost in second pass vs lengths of those retained"
    t2$alternative <- "lengths of lost proteins less than kept"
    p <- t2$`p.value`
    if (p >= 0.05) {
      t2 <- wilcox.test(lengths_lost, lengths_kept, alternative = "greater")
      t2$data.name <- "Lengths of proteins lost in second pass vs lengths of those retained"
      t2$alternative <- "Lengths of lost proteins greater than kept"
    }
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
    if (two_sided_significant != "Y") {
      NA
    } else if (alternative_significant == "Y") {
      glue("first {alt}")
    } else {
      glue("second {alt}")
    }
  }
) |> unlist()

TABLES$pairwise_tests_raw <- pw

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

save(TABLES, glue("{M$outdir}/pass_differences"))
