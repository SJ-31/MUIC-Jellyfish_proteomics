if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
PALETTE <- "ggthemes::colorblind"
PALETTE2 <- "ggthemes::Classic_10_Medium"
library("ggplot2")
TABLES <- list()
GRAPHS <- list()


# File and path setup
open_search_engines <- c("metamorpheusGPTMD", "msfraggerGPTMD", "msfraggerGlyco")
percolator_all <- read_tsv(M$percolator_all)
ENGINES <- percolator_all$engine |> unique()
alignment_types <- c("denovo", "transcriptome", "database", "unmatched_peptide")
standard_search_engines <- c("comet", "identipy", "metamorpheus", "msfragger", "msgf", "tide")


ta <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = ta)
ALIGN_DIR <- glue("{M$outdir}/alignment_metrics")

# Get general alignment metrics
get_general_alignment <- function(path, param, prefix) {
  get_pass <- function(pass) {
    data <- read_tsv(glue("{path}/{pass}/{prefix}_all_wcoverage.tsv")) |>
      mutate(combined_coverage = pcoverage_align) |>
      select(ProteinId, combined_coverage)
    aligned_peptides_path <- glue("{path}/{pass}/aligned_peptides.tsv")
    peptide_map_path <- glue("{path}/{pass}/percolator_peptide_map.tsv")
    per_protein_alignment_metrics_file <- glue("{ALIGN_DIR}/per_protein_{param}_{pass}.tsv")
    if (file.exists(per_protein_alignment_metrics_file)) {
      per_protein_alignment_metrics <- read_tsv(per_protein_alignment_metrics_file)
    } else {
      tracer <- ta$AlignmentTracer(aligned_peptides_path, peptide_map_path)
      per_protein_alignment_metrics <- tracer$run() |>
        as_tibble() |>
        mutate(
          pass = pass, param = param
        ) |>
        inner_join(data, by = join_by(ProteinId))
      write_tsv(
        per_protein_alignment_metrics,
        per_protein_alignment_metrics_file
      )
    }

    combo_alignment_metrics_file <- glue("{ALIGN_DIR}/combos_{param}_{pass}.tsv")
    if (file.exists(combo_alignment_metrics_file)) {
      combo_alignment_metrics <- read_tsv(combo_alignment_metrics_file)
    } else {
      tracer <- ta$AlignmentTracer(aligned_peptides_path, peptide_map_path)
      combo_alignment_metrics <- tracer$run(mode = "combos")
      print(combo_alignment_metrics)
      # if (length(combo_alignment_metrics) > 1) {
      #   combo_alignment_metrics <- combo_alignment_metrics[1]
      # }
      combo_alignment_metrics <- combo_alignment_metrics |>
        as_tibble() |>
        mutate(
          pass = pass, param = param
        ) |>
        inner_join(data, by = join_by(ProteinId))
      write_tsv(
        combo_alignment_metrics,
        combo_alignment_metrics_file
      )
    }

    # Get file for evaluating cost of removing a type of alignment from data
    per_protein_alignment_differences_file <- glue("{ALIGN_DIR}/per_protein_differences_{param}_{pass}.tsv")
    if (file.exists(per_protein_alignment_differences_file)) {
      per_protein_alignment_differences <- read_tsv(per_protein_alignment_differences_file)
    } else {
      tracer <- ta$AlignmentTracer(aligned_peptides_path, peptide_map_path)
      per_protein_alignment_differences <- tracer$run(mode = "differences") |>
        as_tibble() |>
        mutate(
          pass = pass, param = param
        )
      write_tsv(
        per_protein_alignment_differences,
        per_protein_alignment_differences_file
      )
    }
    list(
      metrics = per_protein_alignment_metrics,
      diff = per_protein_alignment_differences
    )
  }
  ran <- lapply(M$passes, get_pass)
  metrics <- bind_rows(
    ran[[1]]$metrics,
    ran[[2]]$metrics
  )
  diff <- bind_rows(
    ran[[1]]$diff,
    ran[[2]]$diff
  )
  return(list(metrics = metrics, diff = diff))
}

all_data <- lapply(seq_along(M$params), \(x) {
  get_general_alignment(M$all_paths[[x]], M$params[[x]], M$prefixes[[x]])
})
all_alignment_metrics <- lapply(all_data, \(x) x$metrics) |> bind_rows()
all_alignment_differences <- lapply(all_data, \(x) x$diff) |> bind_rows()


# ----------------------------------------
# Engine alignments
engine_alignment_metrics <- all_alignment_metrics |> select(
  ProteinId,
  contains(ENGINES), -matches("count|unmatched"), combined_coverage, param, pass
)
TABLES$engine_alignment_metrics <- engine_alignment_metrics


reticulate::source_python(glue("{args$python_source}/plotting.py"))

library("ggpattern")
to_plot <- format_engine_alignment(TABLES$engine_alignment_metrics) |> as_tibble()

GRAPHS$engine_peptide_coverage <- to_plot |>
  ggplot(aes(y = ln_value, x = variable, fill = variable, pattern = type)) +
  geom_boxplot_pattern(size = 1) +
  guides(fill = guide_legend("Engine")) +
  ylab("log coverage (%)") +
  M$default_theme +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  scale_fill_paletteer_d(PALETTE2) +
  scale_pattern_manual(values = c(
    standard = "none", open = "stripe",
    combined = "circle"
  )) +
  scale_x_discrete(
    limits = c("combined", standard_search_engines, open_search_engines)
  ) +
  facet_grid(rows = vars(pass), cols = vars(param))
attr(GRAPHS$engine_peptide_coverage, "width") <- 20


# Identify which engine, if any, is the best-performing

test_helper <- function(cur_param, cur_pass) {
  cov_list <- engine_alignment_metrics %>%
    filter(param == cur_param & pass == cur_pass) |>
    select(-c(ProteinId, pass, param)) |>
    as.list()
  names(cov_list) <- names(cov_list) |> map_chr(\(x) str_replace(x, "_coverage", ""))

  combos <- combn(names(cov_list), 2)

  test_tb <- lapply(
    seq_len(ncol(combos)),
    \(x) {
      greater <- wilcox.test(cov_list[[combos[1, x]]],
        cov_list[[combos[2, x]]],
        alternative = "greater"
      )
      greater$data.name <- glue("{combos[1, x]} x {combos[2, x]}")
      greater$alternative <- glue("{combos[1, x]} greater")
      two_sided <- wilcox.test(
        cov_list[[combos[1, x]]],
        cov_list[[combos[2, x]]]
      )
      two_sided$data.name <- glue("{combos[1, x]} x {combos[2, x]}")
      two_sided$alternative <- glue("two sided")
      bind_rows(htest2tb(greater), htest2tb(two_sided))
    }
  ) %>% bind_rows()
  test_tb |>
    mutate(param = cur_param, pass = cur_pass)
}


all_tests <- lapply(M$params, \(param) {
  lapply(M$passes, \(pass) {
    test_helper(param, pass)
  })
}) |>
  purrr::reduce(\(l, r) {
    bind_rows_list(l, r)
  }) |>
  bind_rows() |>
  get_adjusted_p() |>
  dplyr::select(-c(method, null)) |>
  mutate(pair = paste0(data, " ", param, ",", pass))

sig <- conclude_one_sided(all_tests) |>
  mutate(pass = case_match(pass, "1-First_pass" ~ "1", "2-Second_pass" ~ "2"))

param_str_helper <- function(param, string) {
  nums <- str_extract_all(string, glue("{param}-[12]")) |>
    unlist() |>
    str_remove("[A-Za-z\\-]*")
  glue("{param} ({paste0(nums, collapse = ',')})")
}

TABLES$engine_coverage_pairwise_sig_reduced <- sig |>
  filter(!is.na(conclusion) & !grepl("combined", data)) |>
  group_by(conclusion) |>
  summarise(param = paste0(param, "-", pass, collapse = ";"), count = n()) |>
  mutate(param = map_chr(param, \(pstring) {
    present <- str_split_1(pstring, ";") |>
      map_chr(\(y) str_remove(y, "-.*")) |>
      unique()
    map_chr(present, \(v) param_str_helper(v, pstring)) |> paste0(collapse = ";")
  })) |>
  gt()

TABLES$engine_coverage_pairwise <- gt(all_tests)
TABLES$engine_coverage_pairwise_sig <- sig |> pairwise_conclusion2gt()

# ----------------------------------------
# Evaluate the contribution of each engine onto the protein

engine_tb <- all_alignment_differences |>
  select(ProteinId, pass, param, contains(ENGINES))

engine_longer <- engine_tb |>
  pivot_longer(cols = -c(ProteinId, pass, param))

test_diff_helper <- function(cur_param, cur_pass) {
  current <- engine_tb %>%
    filter(param == cur_param & pass == cur_pass)
  cov_list <- current %>%
    select(-c(ProteinId, pass, param)) |>
    as.list()
  names(cov_list) <- names(cov_list) |> map_chr(\(x) str_replace(x, "_coverage", ""))
  cov_list$total <- current$total_coverage

  tested <- test_all_pairs(cov_list, wilcox.test, two_sided = TRUE) |>
    bind_rows(
      test_all_pairs(cov_list, \(x, y) wilcox.test(x, y, alternative = "less"), alternative_suffix = "less")
    )
  tested |> mutate(param = cur_param, pass = cur_pass)
}

all_diff_tests <- lapply(M$params, \(param) {
  lapply(M$passes, \(pass) {
    test_diff_helper(param, pass) |> rename(data = pair)
  })
}) |>
  purrr::reduce(\(l, r) {
    bind_rows_list(l, r)
  }) |>
  bind_rows() |>
  get_adjusted_p() |>
  dplyr::select(-c(method, null)) |>
  mutate(pair = paste0(data, " ", param, ",", pass))

# "x less than y" means that removing engine x's peptides has a greater impact on coverage compared to y
cov_test_conclusion <- conclude_one_sided(all_diff_tests)
if (nrow(cov_test_conclusion) > 0) {
  cov_test_conclusion <- filter(cov_test_conclusion, !is.na(conclusion))
  TABLES$coverage_impact_conclusion <- gt(cov_test_conclusion)
} else {
  TABLES$no_significant_difference_when_removing_engine_peptides <- 0
}

GRAPHS$engines_removed <- engine_longer |>
  mutate(ln_value = log(value)) |>
  ggplot(aes(y = ln_value, fill = name)) +
  geom_boxplot(size = 1) +
  M$default_theme +
  scale_fill_paletteer_d(PALETTE2) +
  guides(fill = guide_legend(title = "Engine")) +
  ylab("Log coverage (%) if engine's peptides were removed") +
  facet_wrap(~param) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
attr(GRAPHS$engines_removed, "width") <- 20

# ----------------------------------------

# pepmap_calib <- all_pepmaps$Calibrated[1]

# categorize_type <- function(tb) {
#   mutate(tb, type = )
# }


# test_peptides <- function(pepmap_path, param, pass) {
#   pepmap <- read_tsv(pepmap_path)
#   peptide_sets <- get_peptide_sets(pepmap)
#   test_wrapper <- function(var) {
#     sets <- col_from_tb_list(peptide_sets, var) |> discard(\(x) all(is.na(x)))
#     kruskal <- kruskal.test(sets)
#     if (kruskal$p.value > 0.05) {
#       warning(glue("Pairwise tests for peptide variable `{var}` not run, kruskal test not significant"))
#       return(htest2tb(kruskal))
#     }
#     tests <- test_all_pairs(sets, \(x, y) wilcox.test(x, y, alternative = "greater"), "greater") |>
#       bind_rows(test_all_pairs(sets, wilcox.test, two_sided = TRUE)) |>
#       mutate(
#         p_adjust = p.adjust(p_value),
#         significant = ifelse(p_adjust < 0.05, 1, 0)
#       )
#     conclude_one_sided(tests)
#   }

#   # Compare mass, length and ???
#   masses <- col_from_tb_list(peptide_sets, "mass") |> discard(\(x) all(is.na(x)))
#   lengths <- col_from_tb_list(peptide_sets, "length") |> discard(\(x) all(is.na(x)))

#   mass_plot <- gg_numeric_dist(lapply(masses, \(x) log(x)), method = "boxplot") +
#     ylab("log mw") +
#     theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
#     guides(color = "none")

#   length_plot <- gg_numeric_dist(lapply(lengths, \(x) log(x)), method = "boxplot") +
#     ylab("Peptide length") +
#     theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
#     guides(color = "none")

#   mt <- glue("mass_test_result_{param}_{pass}")
#   lt <- glue("length_test_result_{param}_{pass}")
#   TABLES[[mt]] <<- test_wrapper("mass") |> pairwise_conclusion2gt()
#   TABLES[[lt]] <<- test_wrapper("length") |> pairwise_conclusion2gt()
# }

# # Sets of peptides per engine

# # ----------------------------------------
# #' Engine peptide characteristics, per-protein, paired tests

# matched_peptide_map <- {
#   temp <- M$data |>
#     filter(!is.na(MatchedPeptideIds)) |>
#     select(ProteinId, MatchedPeptideIds) |>
#     separate_longer_delim(MatchedPeptideIds, ";")
#   temp
# }

# # Map any matched UP peptides back into DBPs
# pepmap_joined <- left_join(pepmap, matched_peptide_map,
#   by = join_by(x$ProteinId == y$MatchedPeptideIds),
#   relationship = "many-to-many"
# ) |>
#   filter(!is.na(ProteinId.y) | grepl("P", ProteinId)) |>
#   mutate(
#     ProteinId = map2_chr(ProteinId, ProteinId.y, \(x, y) {
#       ifelse(str_detect("P", x) && is.na(y), x, y)
#     })
#   ) |>
#   select(-ProteinId.y) |>
#   filter(!is.na(engine))

# #' Tests engine peptide metrics on a per-protein basis,
# #' rather than by engine groups
# #' @param summarize_fn A function that aggregates the peptide variable `variable`
# #' into a single value for each peptide set per engine
# compare_engine_peptide_metrics <- function(
#     peptide_map,
#     variable = "mass",
#     summarize_fn = \(x) summarize(x, mass = mean(mass))) {
#   combos <- combn(unique(peptide_map$engine), 2)

#   compare_pair <- function(pair) {
#     filtered <- filter(peptide_map, engine %in% pair)
#     cur_engines <- filtered |>
#       group_by(ProteinId) |>
#       filter(length(unique(engine)) == 2) |>
#       group_by(ProteinId, engine) |>
#       summarize_fn() |>
#       ungroup() |>
#       pivot_wider(names_from = engine, values_from = !!as.symbol(variable))
#     x <- cur_engines[[pair[1]]]
#     y <- cur_engines[[pair[2]]]
#     p_string <- glue("{pair[1]} x {pair[2]}")
#     t1 <- wilcox.test(x, y, paired = TRUE) |> htest2tb(data.name = p_string)
#     t2 <- wilcox.test(x, y,
#       paired = TRUE,
#       alternative = "greater"
#     ) |>
#       htest2tb(data.name = p_string, alternative = glue("{pair[1]} greater"))
#     bind_rows(t1, t2)
#   }
#   lapply(seq_len(ncol(combos)), \(i) compare_pair(combos[, i])) |>
#     bind_rows() |>
#     rename(pair = data)
# }


# mass_tests <- compare_engine_peptide_metrics(pepmap_joined) |>
#   get_adjusted_p()
# mass_test_conclusions <- conclude_one_sided(mass_tests)
# length_tests <- compare_engine_peptide_metrics(pepmap_joined,
#   "length",
#   summarize_fn = \(x) summarize(x, length = mean(length))
# ) |>
#   get_adjusted_p()
# length_test_conclusions <- conclude_one_sided(length_tests)


# ----------------------------------------
#' Engine peptide overlaps
#'
#'
all_pepmaps <- lapply(M$all_paths, \(x) {
  c(
    glue("{x}/{M$passes[[1]]}/percolator_peptide_map_all.tsv"),
    glue("{x}/{M$passes[[2]]}/percolator_peptide_map_all.tsv")
  )
}) |> `names<-`(M$params)

# Do this for all runs and parameters
get_overlap_main <- function(peptide_map_path, param, pass) {
  overlap_results <- ta$get_peptide_overlap(peptide_map_path) |> `names<-`(c("overlap_df", "subsets", "peptide_counts"))
  overlap_tb <- overlap_results$overlap_df |> as_tibble()

  if (!dir.exists(glue("{M$outdir}/engine_characteristics"))) {
    dir.create(glue("{M$outdir}/engine_characteristics"))
  }
  outdir <- glue("{M$outdir}/engine_characteristics")
  if (length(overlap_results$subsets) == 0) {
    cat("no subsets",
      file = glue("{outdir}/peptide_subsets_{param}_{pass}.txt")
    )
  } else {
    cat(overlap_results$subsets,
      file = glue("{outdir}/peptide_subsets_{param}_{pass}.txt")
    )
  }
  list(
    counts = as_tibble(overlap_results$peptide_counts) |> mutate(param = param, pass = pass),
    overlap = overlap_tb |> mutate(param = param, pass = pass)
  )
}

overlap_tb <- lapply(names(all_pepmaps), \(p) {
  paths <- all_pepmaps[[p]]
  first <- get_overlap_main(paths[1], p, "First")
  second <- get_overlap_main(paths[2], p, "Second")
  bind_rows_list(first, second)
}) |> purrr::reduce(\(l, r) list(overlap = bind_rows(l$overlap, r$overlap), counts = bind_rows(l$counts, r$counts)))

peptide_counts <- overlap_tb$counts

to_plot <- overlap_tb$overlap |>
  group_by(first, second) |>
  summarise(overlap = mean(overlap)) |>
  ungroup() |>
  mutate(
    pair = mapply(\(x, y) paste0(str_sort(c(x, y)), collapse = "_"), first, second),
    first = gsub("_.*", "", pair),
    second = gsub(".*_", "", pair)
  ) |>
  distinct(pair, .keep_all = TRUE)

pep_counts <- peptide_counts |>
  pivot_longer(cols = -c(param, pass), names_to = "Engine") |>
  group_by(Engine) |>
  summarise(value = mean(value)) |>
  mutate(value = log(value)) |>
  ggplot(aes(x = Engine, y = value, fill = Engine)) +
  geom_bar(stat = "identity") +
  M$default_theme +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("Log n unique peptides") +
  scale_fill_paletteer_d("palettetown::caterpie")

overlap_plot <- to_plot |>
  ggplot(aes(x = first, y = second, fill = overlap)) +
  geom_tile() +
  geom_text(aes(label = round(overlap, 2))) +
  ylab("Engine") +
  xlab("Engine") +
  scale_fill_paletteer_c("ggthemes::Classic Area Green", name = "Overlap coefficient") +
  theme_minimal() +
  M$default_theme +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "inside",
    legend.position.inside = c(.70, .25),
    panel.grid = element_blank(), axis.title.x = element_blank()
  ) +
  scale_x_discrete(position = "top")

GRAPHS$engine_sim_overlap <- cowplot::plot_grid(overlap_plot, pep_counts)
attr(GRAPHS$engine_sim_overlap, "width") <- 18


TABLES$overlap_coefficient <- overlap_tb$overlap
TABLES$peptide_counts <- peptide_counts
save(c(GRAPHS, TABLES), glue("{M$outdir}/engine_characteristics"))
