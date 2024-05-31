SECOND_PASS_ENGINES <- c("identipy", "msgf", "msfragger", "comet")

TABLES <- list()

db <- glue("{PATH}/Databases/decoysWnormal.fasta")
dbSize <- function(fasta_file) {
  seqkitStat(fasta_file)$num_seqs[[1]] %>% as.numeric()
}

complete_db_size <- dbSize(db)
sizes <- lapply(SECOND_PASS_ENGINES, \(x) {
  dbSize(glue("{PATH}/2-Second_pass/BK_databases/{x}_bk_database.fasta"))
}) %>% `names<-`(SECOND_PASS_ENGINES)

getPercolator <- function(pass) {
  dir <- glue("{PATH}/{pass}/Percolator")
  lapply(SECOND_PASS_ENGINES, \(x) {
    read_tsv(glue("{dir}/{x}_percolator_proteins.tsv"))
  }) %>% `names<-`(SECOND_PASS_ENGINES)
}

first_pass <- getPercolator("1-First_pass")
sec_pass <- getPercolator("2-Second_pass")

getStats <- function(engine_list, pass) {
  getRow <- function(lst_subset) {
    tb <- lst_subset[[1]]
    tibble(
      engine = names(lst_subset),
      num_proteins = nrow(tb),
      num_passed_fdr = tb %>% filter(`q-value` < FDR) %>% nrow(),
      num_passed_pep = tb %>% filter(`posterior_error_prob` < FDR) %>% nrow()
    )
  }
  lmap(engine_list, getRow) %>%
    bind_rows() %>%
    mutate(pass = pass)
}

stats <- bind_rows(getStats(first_pass, "first"), getStats(sec_pass, "sec"))
SEQ_MAP <- read_tsv(glue("{PATH}/Databases/seq-header_mappings.tsv"))


ALL_RESULTS <- bind_rows(
  mutate(run$first, pass = "first"),
  mutate(run$second, pass = "sec")
) %>% tibbleDuplicateAt(., "MatchedPeptideIds", ";")
KEPT_PROTEINS <- c(ALL_RESULTS$ProteinId, ALL_RESULTS$MatchedPeptideIds)
ALL_RESULTS <- local({
  expanded <- ALL_RESULTS %>%
    mutate(ProteinId = MatchedPeptideIds) %>%
    filter(!is.na(ProteinId))
  bind_rows(ALL_RESULTS, expanded)
}) %>%
  distinct(ProteinId, .keep_all = TRUE)

# Per-engine analysis function
perEngine <- function(engine) {
  engine_results <- list()
  data <- dplyr::bind_rows(
    tibbleDuplicateAt(first_pass[[engine]], "ProteinId", ",") %>%
      mutate(., pass = "first"),
    tibbleDuplicateAt(sec_pass[[engine]], "ProteinId", ",") %>%
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
  data <- mutate(data, lost_in_sec = !ProteinId %in% sec$ProteinId)

  new_in_sec <- sec %>% filter(!ProteinId %in% first$ProteinId) # Should be empty

  found_in_both <- inner_join(first, sec,
    by = join_by(x$ProteinId == y$ProteinId),
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
  engine_results$pairwise_tests <- pairwiseFromTb(
    found_in_both, to_test, c("less", "greater", "greater"),
    \(x, y, ...) wilcox.test(x, y, paired = TRUE, ...)
  )
  tab <- table(data$from, data$pass) %>% table2Df()
  chi <- chisq.test(tab)
  # Check if the distribution of protein types identified differs between
  # runs. It shouldn't, because the second should only identify proteins in the
  # first
  chi$data.name <- "Frequency of proteins from different sources"
  chi <- htest2Tb(chi)

  kept_data <- data %>% filter(ProteinId %in% KEPT_PROTEINS)
  merged <- inner_join(kept_data, ALL_RESULTS, by = join_by(ProteinId)) %>%
    select(-contains(".x|.y"))


  sources <- c("transcriptome", "denovo", "database")
  engine_results$OR <- lapply(
    sources,
    \(x) {
      table <- table(data$from != x, data$lost_in_sec)
      or <- table %>% oddsRatio()
      upper <- table %>% oddsRatio(CI = TRUE)
      lower <- table %>% oddsRatio(CI = TRUE, side = "lower")
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
  test <- htest2Tb(test)
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
    test <- bind_rows(test, htest2Tb(t2))
  }
  engine_results$htests <- bind_rows(chi, test)
  return(engine_results)
}

TABLES <- purrr::reduce(SECOND_PASS_ENGINES, \(acc, x) {
  mutateBind <- function(col) {
    acc[[col]] <- bind_rows(
      acc[[col]],
      mutate(results[[col]], engine = x)
    )
  }
  results <- perEngine(x)
  acc$pairwise_tests <- mutateBind("pairwise_tests")
  acc$htests <- mutateBind("htests")
  acc$OR <- mutateBind("OR")
  return(acc)
}, .init = list(pairwise_tests = tibble(), htests = tibble(), OR = tibble()))

TABLES$stats <- stats

save(TABLES, glue("{OUTDIR}/figures/pass_differences"))
