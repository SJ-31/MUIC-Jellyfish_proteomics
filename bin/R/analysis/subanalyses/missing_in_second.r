if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
SECOND_PASS_ENGINES <- c("identipy", "msgf", "msfragger", "comet")
library("ggpattern")

# ----------------------------------------
# Compare expect values between passes
get_msgf <- function(dir) {
  to_character <- c("Label")
  to_double <- c("ScanNr", "ExpMass", "CalcMass")
  list.files(dir, pattern = "*pin", full.names = TRUE) |>
    lapply(\(x) {
      read_tsv(x) |> mutate(
        across(all_of(to_character), as.character),
        across(all_of(to_double), as.double)
      )
    }) |>
    bind_rows()
}

JOIN_SUFFIX <- c(".first", ".sec")

SCORE_COLS <- list(
  identipy = "hyperscore",
  comet = "xcorr",
  msfragger = "hyperscore",
  msgf = "RawScore"
)

EXEPCTED_COLS <- list(
  identipy = "expect",
  comet = "e-value",
  msfragger = "expect",
  msgf = "expect"
)

format_pin <- function(tb) {
  filter(tb, Label != -1)
}

format_comet <- function(tb) {
  filter(tb, !grepl("rev_", protein))
}

merge_engines <- function(first, second, engine_name) {
  if (engine_name %in% c("identipy", "msfragger", "msgf")) {
    join_col <- "Peptide"
    if (engine_name == "msfragger") {
      format_fn <- \(x) {
        format_pin(x) |> mutate(
          expect = 10^log10_evalue
        )
      }
    } else if (engine_name == "msgf") {
      format_fn <- \(x) {
        format_pin(x) |> mutate(
          expect = exp(lnEValue)
        )
      }
    } else {
      format_fn <- format_pin
    }
  } else if (engine_name == "comet") {
    format_fn <- format_comet
    join_col <- "modified_peptide"
  }
  rename_fn <- \(x) {
    rename(x,
      e_value = all_of(EXEPCTED_COLS[[engine_name]]),
      psm_score = all_of(SCORE_COLS[[engine_name]])
    )
  }
  f <- format_fn(first) |> rename_fn()
  s <- format_fn(second) |> rename_fn()
  inner_join(f, s, by = join_by(!!as.symbol(join_col)), suffix = JOIN_SUFFIX) |>
    distinct(!!as.symbol(join_col), .keep_all = TRUE)
}

outdir <- glue("{M$outdir}/pass_differences")
results_file <- glue("{outdir}/all_psm_comparisons.tsv")

if (!file.exists(results_file)) {
  psm_comparisons <- tibble()
  for (i in seq_along(M$all_paths)) {
    current_path <- M$all_paths[[i]]
    prefix <- M$prefixes[[i]]
    cur_param <- M$params[[i]]
    passes <- list(First = "1-First_pass", Second = "2-Second_pass")
    engines <- lapply(names(passes), \(x) {
      engine_dir <- glue("{current_path}/{passes[[x]]}/Engines")
      e <- list()
      e$identipy <- read_tsv(glue("{engine_dir}/Identipy/identipy_all_pins.temp"))
      e$comet <- read_tsv(glue("{engine_dir}/Comet/{prefix}_comet.tsv"))
      e$msfragger <- read_tsv(glue("{engine_dir}/MsFragger/fragger_all_pins.temp"))
      e$msgf <- get_msgf(glue("{engine_dir}/MSGF"))
      e
    }) |> `names<-`(names(passes))

    tmp <- lapply(SECOND_PASS_ENGINES, \(x) {
      f <- engines$First[[x]]
      s <- engines$Second[[x]]
      merge_engines(f, s, x) |>
        mutate(engine = x) |>
        select(engine, contains("psm_score"), contains("e_value"))
    }) |>
      bind_rows() |>
      mutate(across(where(is.double), log), param = cur_param)
    psm_comparisons <<- bind_rows(psm_comparisons, tmp)
    write_tsv(psm_comparisons, results_file)
  }
} else {
  psm_comparisons <- read_tsv(results_file)
}

TABLES <- list()
plot <- FALSE
if (plot) {
  reticulate::source_python(glue("{args$python_source}/plotting.py"))
  psm_score_plot <- plotly_psm_comparisons(psm_comparisons, "psm_score")
  plotly_save(psm_score_plot, list(file = glue("{outdir}/plotly_psm_score.png")))
  e_value_plot <- plotly_psm_comparisons(psm_comparisons, "e_value")
  plotly_save(e_value_plot, list(file = glue("{outdir}/plotly_e_value.png")))
}


cols <- c("e_value", "psm_score")
equals <- list(
  engine = c(), var = c(), param = c(), n_equal = c(),
  fr_equal = c()
)
all_tests <- lapply(cols, \(p) {
  test_tb <- lapply(SECOND_PASS_ENGINES, \(x) {
    cur <- psm_comparisons |> filter(engine == x)
    left <- glue("{p}.first")
    right <- glue("{p}.sec")

    param_tb <- lapply(M$params, \(rr)  {
      equals$engine <<- append(equals$engine, x)
      equals$param <<- append(equals$param, rr)
      cur_p <- cur |> filter(param == rr)
      cur_eq <- cur_p |> filter(!!as.symbol(left) == !!as.symbol(right))
      equals$var <<- append(equals$var, p)
      equals$n_equal <<- append(equals$n_equal, nrow(cur_eq))
      equals$fr_equal <<- append(equals$fr_equal, nrow(cur_eq) / nrow(cur_p))

      pair <- glue("First x Second {x}, {p}, {rr}")
      two_sided <- wilcox.test(cur_p[[left]], cur_p[[right]], paired = TRUE) |> htest2tb(data.name = pair)
      greater <- wilcox.test(cur_p[[left]], cur_p[[right]],
        paired = TRUE, alternative = "greater"
      ) |> htest2tb(data.name = pair, alternative = "First greater")
      bind_rows(two_sided, greater) |> mutate(param = rr)
    }) |>
      bind_rows() |>
      mutate(engine = x)

    param_tb
  }) |>
    bind_rows() |>
    mutate(var = p)
}) |>
  bind_rows() |>
  get_adjusted_p() |>
  rename(pair = data)

eq_tb <- as_tibble(equals)



TABLES$psm_evalue_tests <- all_tests |>
  conclude_one_sided() |>
  pairwise_conclusion2gt()

TABLES$psm_evalue_tests_raw <- all_tests



TABLES$equality_graph <- eq_tb |> ggplot(aes(
  y = fr_equal, fill = engine,
  x = engine,
  pattern = var
)) +
  geom_bar_pattern(
    stat = "identity", position = "dodge",
    pattern_density = 0.1,
    pattern_spacing = 0.03
  ) +
  facet_wrap(~param) +
  M$default_theme +
  scale_pattern_manual(values = c(e_value = "none", psm_score = "stripe")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ylab("Proportion") +
  scale_fill_paletteer_d("NineteenEightyR::sonny")

TABLES$equality_graph

save(TABLES, outdir)
