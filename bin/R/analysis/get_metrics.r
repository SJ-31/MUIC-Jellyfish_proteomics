library(tidyverse)
library(ggplot2)
library(ggridges)
library(venn)
library(Peptides)
library(glue)
args <- list(r_source = "./bin/R")
source(glue("{args$r_source}/helpers.r"))

#' Helper function for collecting all the data from a single run for
#' comparison
#'
runData <- function(path) {
  name <- gsub(".*/", "", path)
  results <- list(
    first = read_tsv(glue("{path}/1-First_pass/{name}_all_wcoverage.tsv")),
    sec = read_tsv(glue("{path}/2-Second_pass/{name}_all_wcoverage.tsv"))
  )
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

c_indra <- runData("./results/C_indra")
cov_align <- compareFirstSecL(c_indra, "coverage_alignlen",
                              TRUE, "ProteinId")
num_peps <- compareFirstSecW(c_indra, "num_peps", "ProteinId", TRUE)
test_num_peps <- wilcoxWrapper(num_peps)
test_coverage <- wilcoxWrapper(compareFirstSecW(c_indra, "coverage_alignlen", "ProteinId"))
c_indra_coverage <- passDensityPlot(cov_align, 0.05) + labs(x = "percent coverage")

c_indra_uniques <- passUniques(c_indra)


has_mods <- c_indra$first %>% filter(!is.na(Mods))
unique_mods <- has_mods$Mods %>%
  map(\(x) str_split_1(x, "\\|")) %>%
  unlist() %>%
  map_chr(\(x) gsub(" [1-9]+", "", x)) %>%
  unique() %>%
  discard(\(x) grepl("0$", x))
# Make these into tbl columns
