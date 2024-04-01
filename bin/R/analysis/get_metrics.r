library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)
args <- list(r_source = "./bin/R")
source(glue("{args$r_source}/helpers.r"))
source(glue("{args$r_source}/GO_helpers.r"))


EGGNOG_COLS <- c("EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")

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

splitAndCount <- function(tb, col, pattern) {
  vec <- lapply(tb[[col]], \(x) {
    if (is.na(x)) return(x)
    return(str_split_1(x, pattern))
  }) %>%
    unlist()
  return(getFreqTb(vec, col))
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
  data$lineage <- tb$lineage %>%
    purrr::map_chr(., \(x) parseLineage(x, 3)) %>%
    getFreqTb(., "lineage")
  return(data)
}

#' Get the clade at the specified level from an NCBI lineage string
#' if the level is too specific for the given string, return the most specific clade less than the level
parseLineage <- function(lineage_str, level) {
  if (is.na(lineage_str)) return(NA)
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

#' Helper function for collecting all the data from a single run for
#' comparison
#'
runData <- function(path, remove_dupes) {
  name <- gsub(".*/", "", path)
  results <- list(
    first = read_tsv(glue("{path}/1-First_pass/{name}_all_wcoverage.tsv")),
    sec = read_tsv(glue("{path}/2-Second_pass/{name}_all_wcoverage.tsv"))
  )
  if (!missing(remove_dupes) && remove_dupes) {
    results$first <- results$first %>% distinct(ProteinId, .keep_all = TRUE)
    results$sec <- results$sec %>% distinct(ProteinId, .keep_all = TRUE)
  }
  return(results)
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
                         by = join_by(!!join_on))
    f <- joined[[glue("{compare_col}.x")]]
    s <- joined[[glue("{compare_col}.y")]]
  } else {
    f <- first_sec$first[[compare_col]]
    s <- first_sec$sec[[compare_col]]
  }
  bound <- bind_rows(
    tibble(!!compare_col := f,
           pass = "first"),
    tibble(!!compare_col := s,
           pass = "second"),
  )
  return(bound)
}


notMissing <- function(tb) {
  colSums(!is.na(tb)) %>%
    purrr::map_dbl(., \(x) {
      x / nrow(tb) * 100
    }) %>%
    as.data.frame()
}

#'  Compare first and second passes
#'
#' @description
#' Returns a wide-form tibble containing the data of the first and
#' second passes from "runData". By necessity, must join pass data
compareFirstSecW <- function(first_sec, compare_col, join_on, drop_na) {
  joined <- full_join(first_sec$first,
                      first_sec$sec,
                      by = join_by(!!join_on))
  cols <- paste0(compare_col, c(".x", ".y"))
  chosen <- joined %>%
    select(all_of(c(join_on, cols))) %>%
    rename_with(\(x)
                  ifelse(grepl(".x", x),
                         gsub(".x", ".first", x),
                         gsub(".y", ".sec", x)))
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
      tb <- map %>% as.list() %>% as_tibble()
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
  plot <- ggplot(pass_tb, aes(x = .data[[col]], after_stat(density),
                              color = pass)) +
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
  two_sided <- wilcox.test(first, second, paired = TRUE,
                           na.action = na.omit)
  sec_greater <- wilcox.test(first, second, paired = TRUE,
                             na.action = na.omit,
                             alternative = "less")
  return(list(two_sided = two_sided, sec_greater = sec_greater))
}

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
sv_first <- counts$first$go$GO %>%
  purrr::map_dbl(getSV)
sv_sec <- counts$sec$go$GO %>% purrr::map_dbl(getSV)
sv_compare <- tibble(sv = sv_first, pass = "first") %>%
  bind_rows(tibble(sv = sv_sec, pass = "sec"))

per_protein <- list()
# Counts of go terms per protein were recorded already in combine_all.r

# Max sv per protein (better do this in combine_all







