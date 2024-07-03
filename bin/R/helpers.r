flatten_by <- function(vector, sep, na.rm = TRUE) {
  if (na.rm) {
    vector <- vector |> discard(is.na)
  }
  lapply(vector, \(x) str_split(x, sep)) %>% unlist()
}

clean_peptide <- function(pep) {
  # Leave only residue characters in a peptide sequence
  if (grepl("\\]|[a-z0-9.]|-", pep)) {
    pep <- str_to_upper(pep) %>%
      str_extract_all("[A-Z]") %>%
      unlist() %>%
      paste0(collapse = "")
  }
  return(pep)
}

most_abund <- function(grouped_df) {
  p <- pick(base_peptide)$base_peptide
  up <- unique(p)
  tab <- tabulate(match(p, up))
  abund <- up[tab == max(tab)]
  current <- grouped_df[cur_group_rows(), ]
  choices <- filter(current, base_peptide %in% abund)
  if (nrow(choices) > 1) {
    return(sample_n(choices, 1))
  }
  return(choices)
}

gg_numeric_dist <- function(lst, method = "hist", ...) {
  if (is.null(names(lst))) {
    names(lst) <- LETTERS[seq_along(lst)]
  }
  tb <- purrr::lmap(lst, \(x) tibble(freq = x[[1]], id = names(x))) %>%
    bind_rows()
  if (method == "boxplot") {
    plot <- tb %>% ggplot(aes(x = id, y = freq, color = id)) +
      geom_boxplot(...)
  } else if (method == "hist") {
    plot <- tb %>% ggplot(aes(x = freq, color = id, fill = id)) +
      geom_histogram(...)
  } else if (method == "ridgeline") {
    plot <- tb %>% ggplot(aes(y = id, color = id, fill = id, x = freq)) +
      ggridges::geom_density_ridges(...)
  } else if (method == "freq_poly") {
    plot <- tb %>% ggplot(aes(x = freq, color = id)) +
      geom_freqpoly(...)
  } else {
    warning("No supported method given")
  }
  return(plot)
}

add_prefix <- function(vector, prefix) {
  map_chr(vector, \(x) glue("{prefix}{x}"))
}

#' Given a list A and B (can be named or not) containing the
#' same number of vector elements, return a list C where
#' the corresponding vectors of A and B have been joined
#' @param A list of vectors
#' @param B list of vectors (same number as A)
merge_lists <- function(A, B) {
  if ((!is.null(names(A)) && !is.null(names(B))) && any(names(A) != names(B))) {
    stop("If lists are named, they must have the same elements!")
  } else if (length(A) != length(B)) {
    stop("Lists don't have the same lengths")
  }
  map2(A, B, \(x, y) c(x, y))
}

#' Return a function that for a given string,
#'  substitutes a element in vector `old` with
#' its corresponding element in `new`.
#' `fn` is an additional function to apply regardless of matching
substitute_all <- function(old, new, fn = NULL) {
  fun <- function(x) {
    for (i in seq_along(old)) {
      if (str_detect(x, old[i])) {
        return(gsub(old[i], new[i], x))
      }
    }
    if (!is.null(fn)) {
      return(fn(x))
    }
    return(x)
  }
  fun
}

#' Verify a function's arguments
#' If `arg` fails the predicate, then stop is called, and the name of the
#' function is printed
#' @param arg Function argument
#' @param predicate Predicate function, such as one that checks if `arg`
#' is of a specific class
assertArg <- function(arg, predicate) {
  if (!predicate(arg)) {
    callstack <- sys.calls()
    top_level <- deparse(callstack[[1]])
    fn_start <- str_locate(top_level, "\\(")[, 2] - 1
    caller <- str_sub(top_level, start = 1L, end = fn_start)
    last <- deparse(callstack[[length(callstack)]])
    failed_arg <- str_extract(last, "assertArg\\((.*),", group = 1)
    message <- glue("wrong argument in `{caller}`")
    message <- glue("{message}\n  Incorrect argument `{failed_arg}`")
    stop(message, call. = FALSE)
  }
}

#' Run pairwise tests against elements in `to_test`, using a
#' pre-determined testing function.
#' Note: does not adjust for multiple testing
testAllPairs <- function(to_test, hypothesisTest) {
  assertArg(to_test, \(x) class(x) == "list" && !is.null(names(x)))
  combos <- combn(names(to_test), 2)
  lapply(seq_len(ncol(combos)), \(x) {
    a <- combos[1, x]
    b <- combos[2, x]
    test <- hypothesisTest(
      to_test[[a]],
      to_test[[b]]
    )
    tibble(comparison = glue("{a} x {b}"), p_value = test$p.value, alternative = test$alternative)
  }) %>% bind_rows()
}

table2tb <- function(table, id_col) {
  assertArg(table, \(x) class(x) == "table")
  tb <- as_tibble(table, .name_repair = "unique")
  colnames(tb)[1] <- id_col
  tb
}

table2df <- function(table) {
  assertArg(table, \(x) class(x) == "table")
  convertRow <- function(n) {
    as.list(table[n, ]) %>% as.data.frame()
  }
  purrr::reduce(seq_len(nrow(table)), \(acc, x) {
    bind_rows(acc, convertRow(x))
  }, .init = data.frame()) %>% `rownames<-`(rownames(table))
}



lget <- function(lst, arg, default) {
  find <- lst[[arg]]
  if (is.null(find)) {
    default
  } else {
    find
  }
}

seqkit_stat <- function(filename) {
  call <- system2(
    command = "seqkit", args = c("stat", filename, "-T"),
    stdout = TRUE
  ) %>%
    lapply(., \(x) str_split_1(x, "\t"))
  tb <- as_tibble(setNames(as.list(call[[2]]), call[[1]]))
  if (length(call) > 2) {
    tb <- purrr::reduce(call[3:length(call)], \(x, y) {
      row <- as_tibble(setNames(as.list(y), call[[1]]))
      bind_rows(x, row)
    }, .init = tb)
  }
  return(tb)
}


to <- function(obj, x, val) {
  obj[[x]] <- val
  obj
}

index <- function(vector, a, b = NULL) {
  assertArg(vector, is.atomic)
  if (is.null(b)) {
    return(vector[a, length(vector)])
  }
  return(vector[a:b])
}

reticulate_show_error <- function(expr) {
  captured <- substitute(expr)
  result <- tryCatch(
    expr = eval(captured, envir = parent.frame()),
    error = \(cnd) {
      last_error <- reticulate::py_last_error()
      message("Python error: ", last_error$type, "\n", last_error$value, "\n", last_error$traceback)
    }
  )
  result
}

is_scalar_col <- function(column) {
  all(map_lgl(column, ~ length(.x) == 1))
}

max_peptide_length <- function(peptide_vector) {
  longest <- function(peptide_str) {
    str_split_1(peptide_str, ";") |>
      map_chr(clean_peptide) |>
      nchar() |>
      max()
  }
  map_dbl(peptide_vector, longest)
}

values <- function(vector) {
  copy <- vector
  names(copy) <- NULL
  copy
}

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



DEFAULT_THEME <- theme(
  axis.title.y = element_text(size = 22),
  legend.title = element_text(size = 20, face = "bold")
)

binary_table2tb <- function(table, var) {
  return(tibble(!!var := c(FALSE, TRUE), "FALSE" = table[, 1], "TRUE" = table[, 2]))
}

format_chi_expected <- function(table) {
  bin <- binary_table2tb(table, "var")
  expected <- binary_table2tb(show_expected_chi(table), "var") |>
    mutate(across(is.double, \(x) paste0(" (", round(x, 2), ")")))
  bin %>% mutate(
    `FALSE` = paste0(`FALSE`, expected$`FALSE`),
    `TRUE` = paste0(`TRUE`, expected$`TRUE`)
  )
}


#' Compute odds ratios and chi square tests between a NON mutually exclusive categorical
#' variable `a` and exclusive categorical variable `b`.
#' Note: this also works for exclusive categorical variables `a`
#'
#' Function works around exclusivity requirement of chi square by collapsing contingency
#' table for each level of `a`
#' @param var_a_levels the levels of the non-mutually exclusive categorical variable.
#' `tb` is expected to have columns corresponding to these.
#' @param var_b_col the NON-binary column containing the other categorical var
#' @param tb a semi-binary tb where rows are observations and columns are the values of
#' categorical variables.
#' @param binary If true, levels specified in `var_a_levels` are expected to be binary
#' columns (TRUE/FALSE) indicating when an observation has been flagged with that level
#'    Otherwise, you must specify `var_a_col`, where its entries are the given level for #'   `var_a`.
chisqNME <- function(
    tb, var_a_levels, var_b_col, var_a, var_b, binary = TRUE,
    var_a_col = NULL) {
  gt <- list()
  tables <- list()
  var_b_levels <- unique(tb[[var_b_col]])
  chi <- tibble()
  table_list <- list()

  for (a_level in var_a_levels) {
    table_list[[a_level]] <- list()
    for (b_level in var_b_levels) {
      if (binary) {
        t <- table(tb[[var_b_col]] == b_level, tb[[a_level]])
      } else {
        t <- table(tb[[var_b_col]] == b_level, tb[[var_a_col]] == a_level)
      }
      table_list[[a_level]][[b_level]] <- t
      tst <- tryCatch(
        expr = chisq.test(t),
        error = \(cnd) {
          print(glue("a_level: {a_level}, b_level: {b_level}"))
          print(t)
          stop("Chi square error")
        }
      )
      row <- tibble(!!var_a := a_level, !!var_b := b_level, p_value = tst$p.value)
      chi <- bind_rows(chi, row)
    }
  }
  chi <- chi %>% mutate(
    p_adjust = p.adjust(p_value, "holm"),
    is_significant = ifelse(p_value < 0.05, TRUE, FALSE)
  )

  # Compute effect size for significant hits, using the odds ratio
  # This represents how much more likely the engine is to identify
  # peptides from a protein of the given category
  # Confidence interval is 95%
  odds_ratios <- chi %>%
    apply(1, \(x) {
      row <- tibble(OR = NA, OR_upper_ci = NA, OR_lower_ci = NA)
      if (!as.logical(x["is_significant"])) {
        return(row)
      }
      e <- x[var_a]
      c <- x[var_b]
      table <- table_list[[e]][[c]]
      try({
        odds_ratio <- get_odds_ratio(table)
        odds_ratio_ci_u <- get_odds_ratio(table, TRUE, "upper")
        odds_ratio_ci_l <- get_odds_ratio(table, TRUE, "lower")
        row <- tibble(
          OR = odds_ratio,
          OR_upper_ci = odds_ratio_ci_u,
          OR_lower_ci = odds_ratio_ci_l
        )
      })
      return(row)
    }) %>%
    bind_rows()

  chi_f <- bind_cols(chi, odds_ratios)
  tables$chi <- chi_f
  gt$chi <- chi_f %>%
    mutate(across(is.double, \(x) round(x, 4))) %>%
    mutate(
      OR = paste0(OR, " [", OR_lower_ci, ", ", OR_upper_ci, "]"),
      OR = map_chr(OR, \(x) ifelse(grepl("NA|NaN", x), "NA", x))
    ) %>%
    select(-matches("OR_")) %>%
    rename("Odds ratio, 95% CI [lower, upper]" = OR) %>%
    gt() %>%
    tab_header(
      title = glue("Association between {var_a} and {var_b}"),
      subtitle = "Tested using chi-square"
    )
  # Record all results
  var_a_table_tb <- tibble(!!var_b := var_b_levels)
  remove_cat <- FALSE
  for (e in names(table_list)) {
    e_list <- table_list[[e]]
    temp_lst <- list()
    remove_col_label <- FALSE
    for (ct in names(e_list)) {
      tib <- format_chi_expected(e_list[[ct]]) |> rename(!!var_b := var)
      if (remove_cat) {
        tib <- dplyr::select(tib, -!!var_b)
      }
      tib <- gt(tib)
      if (!remove_col_label) {
        remove_col_label <- TRUE
      } else {
        tib <- tib %>% tab_options(column_labels.hidden = TRUE)
      }
      temp_lst <- c(temp_lst, as_raw_html(tib))
    }
    temp_tb <- tibble(!!e := temp_lst)
    var_a_table_tb <- bind_cols(var_a_table_tb, temp_tb)
    remove_cat <- TRUE
  }
  tables$contingency <- var_a_table_tb
  gt$contingency <- var_a_table_tb %>%
    gt() %>%
    fmt_markdown(columns = everything()) %>%
    tab_header(
      title = md(glue("**Contingency tables used in chi-square analysis of {var_a} vs {var_b}**")),
    )
  return(list(gt = gt, tb = tables))
}

#' Alternative method of grouping protein identifications by their shared unique peptides
group_by_unique_peptides <- function(tb) {
  if ("GroupUP" %in% colnames(tb)) {
    print("Has Group already, ignoring")
    return(tb)
  }
  split_with_na <- function(str) {
    if (is.na(str)) {
      return("NA")
    }
    str_split_1(str, ";")
  }
  split_tb <- tb |> mutate(
    unique_peptides =
      lapply(unique_peptides, split_with_na)
  )
  new_groups <- split_tb %>%
    group_by(unique_peptides) |>
    nest()
  new_groups$size <- map_dbl(new_groups$data, nrow)
  new_groups$GroupUP <- paste0("G", seq_len(nrow(new_groups)))
  group_map <- new_groups |>
    unnest(cols = c(data)) |>
    ungroup() |>
    select(ProteinId, GroupUP)
  tb |>
    inner_join(group_map) |>
    relocate(GroupUP, .after = Group)
}

#' Replace discrete labels in a ggplot
#'
#' @description
#' @param mapping a named vector where names are the new labels and values are the old
#' @value a function that takes in breaks and converts the old labels to the new ones
#' Use this as "scale_x_discrete(labels = get_label_replacement(YOUR_LABEL_MAPPING))"
get_label_replacement <- function(mapping) {
  \(breaks) {
    reversed <- as.list(setNames(names(mapping), mapping))
    map_chr(breaks, \(x) {
      if (x %in% names(reversed)) {
        reversed[[x]]
      } else {
        x
      }
    })
  }
}

#' Automatically determine which variable was greater/lesser from a tibble of the
#' one-sided hypothesis tests
#'
#' @description
#' The assumption is that the two sided test for this pair is significant
get_one_sided_conclusion <- function(alternative, pair, pair_sep, is_significant) {
  if (alternative == "two sided") {
    return(NA)
  }
  splits <- str_split_1(pair, pair_sep) |> map_chr(str_trim)
  first <- splits[1]
  second <- splits[2]
  direction <- case_when(
    str_detect(alternative, "greater") ~ "greater",
    str_detect(alternative, "less") ~ "less",
    .default = NA
  )
  if (is.na(direction)) {
    stop("The alternative must be formulated as \"`x` greater/less!\"")
  }
  alternative_var <- str_replace(alternative, direction, "") |> str_trim()
  if (!alternative_var %in% splits) {
    stop("No recognizable variable for the alternative specified! It must be one of the pair")
  }
  case_when(
    is_significant == 1 && alternative_var == first ~ glue("{alternative_var} {direction} than {second}"),
    is_significant == 1 && alternative_var == second ~ glue("{alternative_var} {direction} than {first}"),
    is_significant == 0 && alternative_var == first ~ glue("{second} {direction} than {first}"),
    is_significant == 0 && alternative_var == second ~ glue("{first} {direction} than {first}"),
  )
}

#' Conclude the results of a series of one and two-sided tests
#'
#' @description
#' Automatically determine which variable in a one-sided test had the statistically
#' significant direction
#' @value a tibble containing pairs that were significant, and a `conclusion` column
#' indicating how to interpret the one-sided test result
conclude_one_sided <- function(htest_tb, pair_sep = "x") {
  unique_pairs <- htest_tb$pair |> unique()
  two_sided_was_significant <- lapply(unique_pairs, \(x) {
    row <- htest_tb |> filter(pair == x & alternative == "two sided")
    if (row$significant == 1) {
      return(x)
    }
    NA
  }) |>
    unlist() |>
    discard(is.na)
  tb <- htest_tb |> filter(pair %in% two_sided_was_significant)
  winning_vars <- pmap(
    list(tb$alternative, tb$pair, rep(pair_sep, length(tb$pair)), tb$significant),
    get_one_sided_conclusion
  ) |> unlist()
  tb |>
    mutate(conclusion = winning_vars) |>
    filter(!is.na(conclusion))
}
