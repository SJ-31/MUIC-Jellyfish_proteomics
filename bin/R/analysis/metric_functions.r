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
runData <- function(prefix, remove_dupes, path) {
  results <- list(
    first = read_tsv(glue("{path}/1-First_pass/{prefix}_all_wcoverage.tsv")),
    sec = read_tsv(glue("{path}/2-Second_pass/{prefix}_all_wcoverage.tsv"))
  )
  dl_f <- getDeeploc(
    glue("{path}/1-First_pass/Deeploc/deeploc_results.csv"),
    glue("{path}/1-First_pass/Combined/unified_groups.tsv")
  )
  dl_s <- getDeeploc(
    glue("{path}/1-First_pass/Deeploc/deeploc_results.csv"),
    glue("{path}/1-First_pass/Combined/unified_groups.tsv")
  )
  results$first <- bind_rows(results$first, dl_f)
  results$sec <- bind_rows(results$sec, dl_s)
  if (!missing(remove_dupes) && remove_dupes) {
    results$first <- results$first %>% distinct(ProteinId, .keep_all = TRUE)
    results$sec <- results$sec %>% distinct(ProteinId, .keep_all = TRUE)
  }
  return(results)
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
  for (col in EGGNOG_COLS) {
    data[[col]] <- splitAndCount(tb, col, ";|,")
  }
  data$PANTHER <- splitAndCount(tb, "PANTHER", ";")
  data$go <- goVector(tb, go_column = "GO_IDs") %>% getFreqTb(., "GO")
  data$matched_peptides <- splitAndCount(tb, "MatchedPeptideIds", ";",
                                         func = removeDigits, unique_only = TRUE
  )
  data$lineage <- tb$lineage %>%
    purrr::map_chr(., \(x) parseLineage(x, 3)) %>%
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
  has_mods <- dplyr::filter(run_df, !is.na((Mods)))
  lst <- list()
  df <- has_mods$Mods %>%
    # Store modification info in a map of mod name -> count
    lapply(., \(m) {
      map <- hash::hash()
      str_split_1(m, "\\|") %>% lapply(., \(x) {
        split <- str_split_1(x, " ")
        map[split[1]] <- split[2]
      })
      tb <- map %>%
        as.list() %>%
        as_tibble()
      return(tb)
    }) %>%
    bind_rows() %>%
    mutate(ProteinId = has_mods$ProteinId) %>%
    relocate(ProteinId, .before = everything()) %>%
    mutate(across(!contains("ProteinId"), as.numeric)) %>%
    column_to_rownames(var = "ProteinId")
  mod_percent <- colSums(df, na.rm = TRUE) / sum(df, na.rm = TRUE) * 100
  lst$count_df <- df
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

#' Perform Wilcoxon signed-rank test between passes
#'
#' @description
#' Wrapper function for making tests more convenient, needs wide-format
#' data
wilcoxWrapper <- function(first_sec) {
  col <- first_sec %>%
    colnames() %>%
    keep(\(x) grepl(".first", x)) %>%
    gsub(pattern = ".first", replacement = "", .)
  first <- first_sec[[glue("{col}.first")]]
  second <- first_sec[[glue("{col}.sec")]]
  two_sided <- wilcox.test(first, second,
                           paired = TRUE,
                           na.action = na.omit
  )
  sec_greater <- wilcox.test(first, second,
                             paired = TRUE,
                             na.action = na.omit,
                             alternative = "less"
  )
  return(list(two_sided = two_sided, sec_greater = sec_greater))
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


replaceNaAll <- function(df, value = 0) {
  df %>% mutate_all(~replace(., is.na(.), value))
}
