JOIN_SUFFIX <- c(".first", ".sec")


#' Helper for collecting deep loc results and combining it with additional metadata
getDeeploc <- function(deeploc_path, unmatched_path) {
  dl <- read_csv(deeploc_path)
  um <- read_tsv(unmatched_path)
  merged <- inner_join(um, dl,
    by = join_by(x$ProteinId == y$Protein_ID)
  ) %>%
    mutate(
      inferred_by = "deeploc",
      localization = Localizations,
      category = Signals,
    ) %>%
    select(all_of(c(colnames(um), "inferred_by", "localization", "category")))
  return(merged)
}


#' Helper function for collecting all the data from a single run for
#' comparison
#'
runData <- function(prefix, path, remove_dupes = TRUE, which = "both") {
  passData <- function(pass) {
    tb <- read_tsv(glue("{path}/{pass}/{prefix}_all_wcoverage.tsv")) %>% filter(q_adjust < FDR)
    deeploc <- getDeeploc(
      glue("{path}/{pass}/Deeploc/deeploc_results.csv"),
      glue("{path}/{pass}/Combined/intersected_searches.tsv")
    )
    tb <- bind_rows(tb, deeploc)
    if (remove_dupes) {
      tb <- tb %>% distinct(ProteinId, .keep_all = TRUE)
    }
    return(tb)
  }
  if (which == "both") {
    return(list(first = passData("1-First_pass"), second = passData("2-Second_pass")))
  } else if (which == "first" || which == "1-First_pass") {
    return(passData("1-First_pass"))
  } else {
    return(passData("2-Second_pass"))
  }
}

alignmentData <- function(path, which = "combine") {
  helper <- function(pass) {
    lst <- list(
      replacements = read_tsv(glue("{path}/{pass}/all_replacements.tsv")),
      metrics = read_tsv(glue("{path}/{pass}/alignment_metrics.tsv")),
      peptides = read_tsv(glue("{path}/{pass}/aligned_peptides.tsv"))
    )
    return(lst)
  }
  passes <- c("1-First_pass", "2-Second_pass")
  if (which == "combine") {
    names <- c("first", "sec")
    result <- list(
      replacements = tibble(),
      metrics = tibble(),
      peptides = tibble()
    )
    for (i in c(1, 2)) {
      lst <- helper(passes[i]) |> lapply(, \(x) mutate(x, pass = names[i]))
      for (i in names(result)) {
        result[[i]] <- bind_rows(result[[i]], lst[[i]])
      }
    }
    return(result)
  } else if (which == passes[1]) {
    return(helper(passes[1]))
  } else {
    return(helper(passes[2]))
  }
}


notMissing <- function(tb) {
  colSums(!is.na(tb)) %>%
    purrr::map_dbl(., \(x) {
      x / nrow(tb) * 100
    }) %>%
    as.data.frame()
}


#' Return a tibble of frequencies from the vector `vec`, renaming the first column
#' with `name`
#'
getFreqTb <- function(vec, name) {
  tb <- vec %>%
    discard(is.na) %>%
    table() %>%
    sort(decreasing = TRUE) %>%
    as_tibble()
  if (nrow(tb) > 0) {
    return(rename(tb, !!name := "."))
  }
  return(NA)
}

#' Get the clade at the specified level from an NCBI lineage string
#' if the level is too specific for the given string, return the most specific clade less than the level
parseLineage <- function(lineage_str, level) {
  if (is.na(lineage_str)) {
    return(NA)
  }
  split <- str_split_1(lineage_str, ";")
  len <- length(split)
  if (len >= level) {
    return(split[level])
  }
  while (len < level) {
    level <- level - 1
  }
  return(split[level])
}

#' For a specified `col` of `tb`, where elements are concatenated by a separator
#' specified with `pattern`, split the elements then return a tibble
#' of their frequencies. Apply apply `func` to elements if supplied
splitAndCount <- function(tb, col, pattern, unique_only, func = NULL) {
  if (!missing(unique_only) && unique_only) {
    uniq <- TRUE
  } else {
    uniq <- FALSE
  }
  vec <- lapply(tb[[col]], \(x) {
    if (is.na(x)) {
      return(x)
    }
    if (is.logical(x)) {
      browser()
    }
    split <- str_split_1(x, pattern)
    if (uniq) {
      split <- unique(split)
    }
    if (!is.null(func)) {
      split <- lapply(split, func) %>% unlist()
    }
    return(split)
  }) %>%
    unlist()
  return(getFreqTb(vec, col))
}


removeDigits <- function(x) {
  return(gsub("[0-9]+", "", x))
}

#' Count the frequencies of each element in the columns
#'
getCounts <- function(tb) {
  data <- list()
  data$category <- getFreqTb(tb$category, "category")
  data$organism <- getFreqTb(tb$organism, "organism")
  # TODO: Restore this once you've figured out the error with
  # the columns being logicals
  # for (col in EGGNOG_COLS) {
  # if (col in colnames(tb)) {
  #   data[[col]] <- splitAndCount(tb, col, ";|,")
  # }
  # }
  data$PANTHER <- splitAndCount(tb, "PANTHER", ";")
  data$go <- goVector(tb, go_column = "GO_IDs") %>% getFreqTb(., "GO")
  data$matched_peptides <- splitAndCount(tb, "MatchedPeptideIds", ";",
    func = removeDigits, unique_only = TRUE
  )
  data$lineage <- tb$lineage %>%
    purrr::map_chr(., \(x) parseLineage(x, 5)) %>%
    getFreqTb(., "lineage")
  return(data)
}

#'  Compare first and second passes
#'
#' @description
#' Returns a long-form tibble containing the data of the first and
#' second passes from "runData"
compareFirstSecL <- function(first_sec, compare_col, compare_found, join_on) {
  if (!missing(compare_found) && compare_found) {
    if (missing(join_on)) {
      stop('Missing "join_on" argument')
    }
    joined <- inner_join(first_sec$first,
      first_sec$sec,
      by = join_by(!!join_on)
    )
    f <- joined[[glue("{compare_col}.x")]]
    s <- joined[[glue("{compare_col}.y")]]
  } else {
    f <- first_sec$first[[compare_col]]
    s <- first_sec$sec[[compare_col]]
  }
  bound <- bind_rows(
    tibble(!!compare_col := f,
      pass = "first"
    ),
    tibble(!!compare_col := s,
      pass = "second"
    ),
  )
  return(bound)
}

#'  Compare first and second passes
#'
#' @description
#' Returns a wide-form tibble containing the data of the first and
#' second passes from "runData". By necessity, must join pass data
compareFirstSecW <- function(first_sec, compare_col, join_on, drop_na) {
  joined <- full_join(first_sec$first,
    first_sec$sec,
    by = join_by(!!join_on),
    suffix = JOIN_SUFFIX
  )
  cols <- paste0(compare_col, JOIN_SUFFIX)
  chosen <- joined %>%
    select(all_of(c(join_on, cols)))
  if (!missing(drop_na) && drop_na) {
    chosen <- filter(chosen, !if_any(everything(), is.na))
  }
  return(chosen)
}

#' Filter proteins unique to each pass
#'
passUniques <- function(first_sec) {
  sec <- first_sec$sec %>%
    filter(!ProteinId %in% first_sec$first$ProteinId)
  first <- first_sec$first %>%
    filter(!ProteinId %in% first_sec$sec$ProteinId)
  return(list(first = first, sec = sec))
}

#' Get protein modification metrics
#'
#' @description
#' Parses the ptm strings in the run tibble to count the occurence of each
#' ptm. Returns...
#' 1. dataframe of protein ids as the index and ptms as columns
#' 2. Proportion of each ptm in the given sample
modMetrics <- function(run_df) {
  getMetricTb <- function(mod_str) {
    map <- hash::hash()
    str_split_1(mod_str, "\\|") |> lapply(\(x) {
      split <- str_split_1(x, " ")
      map[split[1]] <- split[2]
    })
    tb <- map %>%
      as.list() %>%
      as_tibble()
    return(tb)
  }
  has_mods <- dplyr::filter(run_df, !is.na((Mods)))
  lst <- list()
  df <- has_mods$Mods %>%
    # Store modification info in a map of mod name -> count
    lapply(., getMetricTb) %>%
    bind_rows() %>%
    mutate(ProteinId = has_mods$ProteinId) %>%
    relocate(ProteinId, .before = everything()) %>%
    mutate(across(!contains("ProteinId"), as.numeric)) %>%
    column_to_rownames(var = "ProteinId")
  mod_percent <- colSums(df, na.rm = TRUE) / sum(df, na.rm = TRUE) * 100
  lst$count_df <- df |> replaceNaAll(0)
  lst$percentages <- mod_percent
  return(lst)
}

#' Plot the densities of a select column from each pass
#'
#' @arg pass_tb a long-form tibble, ideally the output of compareFirstSecond
passDensityPlot <- function(pass_tb, bw) {
  col <- pass_tb %>%
    colnames() %>%
    discard(\(x) x == "pass")
  plot <- ggplot(pass_tb, aes(
    x = .data[[col]], after_stat(density),
    color = pass
  )) +
    geom_freqpoly(binwidth = bw)
  return(plot)
}

#' Perform Wilcoxon signed-rank test between passes: x is first, y is
#' second
#'
#' @description
#' Wrapper function for making tests more convenient, needs wide-format
#' data
wilcoxWrapper <- function(first_sec, paired, metric) {
  col <- first_sec %>%
    colnames() %>%
    keep(\(x) grepl(".first", x)) %>%
    gsub(pattern = ".first", replacement = "", .)
  first <- first_sec[[glue("{col}.first")]]
  second <- first_sec[[glue("{col}.sec")]]
  tests <- list()
  hypotheses <- c("less", "greater")
  for (h in hypotheses) {
    tests[[h]] <- wilcox.test(first, second,
      paired = paired,
      na.action = na.omit,
      alternative = h
    )
  }
  tb <- tibble(
    alternative = hypotheses,
    metric = rep(metric, 3),
    p_value = map_dbl(tests, \(x) x$p.value),
    statistic = map_dbl(tests, \(x) x$statistic)
  ) %>%
    mutate(reject_null = map_chr(
      p_value,
      \(x) ifelse(x < 0.05, "Y", "N")
    ))
  return(tb)
}


#' Wrapper function for performing pairwise tests
pairwiseFromTb <- function(
    tb, compare_cols, one_sided_alts, testFun, alpha = 0.05,
    suffixes = c(".first", ".sec")) {
  if (length(compare_cols) != length(one_sided_alts)) {
    stop("Must provide alternative hypothesis for each pair of columns!")
  }
  map2(
    compare_cols,
    one_sided_alts,
    \(x, y) {
      f <- tb[[glue("{x}{suffixes[1]}")]]
      s <- tb[[glue("{x}{suffixes[2]}")]]
      ts_test <- testFun(f, s)
      os_test <- testFun(f, s, alternative = y)
      mean_diff <- mean(s, na.rm = TRUE) - mean(f, na.rm = TRUE)
      tibble(
        metric = x, two_sided_p_value = ts_test$`p.value`,
        one_sided_p_value = os_test$`p.value`,
        alternative = glue("{suffixes[1]} {y}"),
        mean_diff = mean_diff
      ) %>% mutate(
        two_sided_significant = map_chr(two_sided_p_value, \(x) {
          ifelse(x < alpha, "Y", "N")
        }),
        alternative_significant = map2_chr(two_sided_p_value, one_sided_p_value, \(x, y)  {
          case_when(
            x < alpha && y < alpha ~ "Y",
            x < alpha && y > alpha ~ "N",
            x > alpha ~ "N"
          )
        })
      )
    }
  ) %>%
    bind_rows()
}


#' Combine lfq results from flashlfq, maxlfq and directlfq, which
#' returns log2 intensity across the different tools
#'
#' @description
#' @param target the metric to average, either "mean" or "median"
mergeLfq <- function(tb, target) {
  lfq_cols <- paste0(c("directlfq", "maxlfq", "flashlfq"), "_", target)
  quant <- tb %>% dplyr::select(ProteinId, contains(lfq_cols))
  is.na(quant) <- do.call(cbind, lapply(quant, is.infinite))
  quant <- mutate(quant,
    log_intensity = pmap(
      list(maxlfq_mean, directlfq_mean, flashlfq_mean),
      \(x, y, z) mean(c(x, y, z), na.rm = TRUE)
    ) %>% unlist()
  )
  return(quant)
}

#' Combines lfq results more flexibly, considering each lfq column
#' rather than their means/medians
#'
#' @description
#' @param combineFun A function that receives a tibble/df as its only
#' argument
mergeLfqCustom <- function(tb, combineFun) {
  quant <- tb %>%
    select(matches("directlfq|flashlfq|maxlfq")) %>%
    select(-matches("mean|median"))
  is.na(quant) <- do.call(cbind, lapply(quant, is.infinite))
  return(combineFun(quant))
}


replaceNaAll <- function(df, value = 0) {
  df %>% mutate_all(~ replace(., is.na(.), value))
}

minMaxScaler <- function(vec) {
  max <- max(vec, na.rm = TRUE)
  min <- min(vec, na.rm = TRUE)
  return(purrr::map_dbl(vec, \(x) (x - min) / (max - min)))
}

#' Returns
#' 1. a list mapping KEGG pathways to the proteins of `tb` that participate
#' in them
#' 2. A tibble mapping ProteinIds->KEGG Pathways, used to generate 1
#' Optionally filter pathways that have fewer than the `minimum` proteins
groupPathways <- function(tb, minimum = 20) {
  id2pathway <- tb %>%
    filter(!is.na(KEGG_Pathway)) %>%
    dplyr::select(ProteinId, KEGG_Pathway) %>%
    tibbleDuplicateAt("KEGG_Pathway", "[;,]")
  pathway_lists <- id2pathway %>%
    group_by(KEGG_Pathway) %>%
    nest() %>%
    apply(1, \(x) {
      lst <- list()
      lst[[x$KEGG_Pathway]] <- x$data$ProteinId
      return(lst)
    }) %>%
    do.call(c, .) %>%
    discard(\(x) length(x) < minimum)
  return(list(map = id2pathway, grouped = pathway_lists))
}

#' Split and a ProteinGroupId string by the ";", optionally remove the
#' numbers and leave unique groups
#'
splitGroupStr <- function(group_str, remove_nums = FALSE, unique = FALSE) {
  if (remove_nums) {
    group_str <- gsub("[0-9]+", "", group_str)
  }
  split <- str_split_1(group_str, ";")
  if (unique) split <- base::unique(split)
  return(split)
}

#' For a nxm contingency table, return a table depicting the expected values used
#' by the chi square test of independence
showExpectedChi <- function(table) {
  mat <- matrix(table, nrow = nrow(table), ncol = ncol(table))
  s <- sum(mat)
  ct <- colSums(mat)
  rt <- rowSums(mat)
  for (i in seq_len(ncol(mat))) {
    for (j in seq_len(nrow(mat))) {
      mat[j, i] <- (ct[i] * rt[j]) / s
    }
  }
  table <- as.table(mat) %>%
    `colnames<-`(colnames(table)) %>%
    `rownames<-`(rownames(table))
  return(table)
}


formatEngineContingency <- function(table, category, w_expected = FALSE) {
  formatted <- table %>%
    as_tibble(.name_repair = "unique") %>%
    pivot_wider(names_from = `...2`, values_from = n)
  if (!missing(category)) {
    formatted <- dplyr::mutate(formatted, `...1` = c(glue("not {category}"), glue("is {category}")), .before = dplyr::everything())
  } else {
    formatted <- dplyr::mutate(formatted, `...1` = c("FALSE", "TRUE"), .before = dplyr::everything())
  }
  formatted <- dplyr::rename(formatted,
    was_hit = "TRUE",
    not_hit = "FALSE", "category" = `...1`
  )
  if (!w_expected) {
    return(formatted)
  }
  expected <- formatEngineContingency(
    showExpectedChi(table),
    FALSE
  ) %>%
    mutate(across(is.double, \(x) paste0(" (", round(x, 2), ")")))
  return(formatted %>% mutate(
    not_hit = paste0(not_hit, expected$not_hit),
    was_hit = paste0(was_hit, expected$was_hit)
  ))
}

#' Compute the odds ratio from a 2x2 contingency table, or the
#' upper and lower 95% CIs of that ratio
#' Table should be of the standard epidemiologic form
#'                Present Absent
#' Exposed          a       b
#' Not exposed      c       d
#' The calculated ratio (OR) is then the odds of `Present` is OR times as
#' high in `Exposed` compared to `Not Exposed`
oddsRatio <- function(ctable, CI = FALSE, side = "upper") {
  a <- ctable[1, 1]
  d <- ctable[2, 2]
  b <- ctable[1, 2]
  c <- ctable[2, 1]
  odds_ratio <- (a * d) / (c * b)
  if (!CI) {
    return(odds_ratio)
  }
  if (side == "upper") {
    x <- log(odds_ratio) + 1.96 * sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  } else {
    x <- log(odds_ratio) - 1.96 * sqrt(1 / a + 1 / b + 1 / c + 1 / d)
  }
  return(exp(1)^x)
}


#' Get a list of vectors in specific groups
#'
#' @description
#' Given a vector of unique values `v` from a column `col_from`, this function
#' groups the tibble `tb` by each value, and extracts the all values from a specified
#' column `target_col`. The result is a named list grouping the values of `target_col`
#' by `v`
groupListFromTb <- function(tb, v, col_from, target_col) {
  list <- lapply(v, \(x) dplyr::filter(tb, !!as.symbol(col_from) == x) %>%
    purrr::pluck(target_col) %>%
    discard(is.na)) %>% `names<-`(v)
  return(list)
}

avgStdevs <- function(tb, cols, stat) {
  result <- vector()
  for (col in cols) {
    result <- c(result, stat(tb[[col]]))
  }
  return(result)
}


htest2Tb <- function(test) {
  row <- tibble(
    null = test$`null.value`,
    alternative = test$alternative,
    method = test$method,
    data = test$`data.name`,
    statistic = test$statistic,
    p_value = test$`p.value`,
    estimate = test$estimate,
  )
  estimate_name <- names(row$estimate)
  if (!is.null(estimate_name)) {
    row$estimate_name <- estimate_name
  }
  if (!is.null(test$conf.int)) {
    row$lower_ci <- test$conf.int[1]
    row$upper_ci <- test$conf.int[2]
    level <- attr(test$conf.int, "conf.level")
    if (!is.null(level)) {
      row$conf_level <- level
    }
  }
  row
}

tbTranspose <- function(tb) {
  df <- tb %>%
    t() %>%
    as.data.frame()
  colnames(df) <- df[1, ]
  df[-1, ] %>% mutate(across(dplyr::everything(), as.numeric))
}
