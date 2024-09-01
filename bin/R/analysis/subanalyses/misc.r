library("paletteer")
library("ggVennDiagram")
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
TABLES <- list()
GRAPHS <- list()
# --------------------------------------------------------
# Investigating trends in missing quantification
current <- M$data |> inner_join(M$lfq, by = join_by(ProteinId))
cq <- c("directlfq", "maxlfq")
missing_quant_tests <- list()
for (q in cq) {
  noq <- current %>% filter(is.na(!!as.symbol(glue("{cq}_mean"))))
  hasq <- current %>% filter(!ProteinId %in% noq$ProteinId)
  missing_quant_tests[[q]] <- wilcox.test(noq$num_peps, hasq$num_peps,
    alternative = "l"
  )
}
capture.output(missing_quant_tests, file = glue("{M$outdir}/missing_quantification_tests.txt"))
rm(noq)
rm(hasq)
# --------------------------------------------------------
# Get Count of fragments

fragment_names <- c("fragment", "partial")
fragment_regex <- paste0(fragment_names, collapse = "|")

data <- read_tsv(M$data_w_cat_path) |> mutate(
  is_fragment =
    as.double(str_detect(str_to_lower(entry_name), fragment_regex)),
  is_denovo = str_detect(ProteinId, "D")
)

fragments <- data |> filter(is_fragment == 1)

grouped <- data |>
  group_by(GroupUP) |>
  summarize(
    assigned_COG = paste0(unique(assigned_COG), collapse = ";"),
    entry_name = paste0(entry_name, collapse = ";"),
    is_fragment = sum(is_fragment),
    is_denovo = sum(is_denovo),
    size = n()
  ) |>
  mutate(entry_name = map_chr(entry_name, split_unique_join)) |>
  arrange(desc(size))


TABLES$grouped_cog_sizes <- grouped

TABLES$fragment_info <- glue("
Total number of groups: {nrow(grouped)}
Number of protein fragments: {nrow(fragments)}
Number of groups with no fragments: {(grouped$is_fragment == 0) |> sum()}
Number of groups consisting only of fragments: {filter(grouped, is_fragment == size) |> nrow()}
Number of groups containing fragments: {length(unique(fragments$GroupUP))}
")



# Just to check if group assignments are correct
cog_evidence <- read_tsv(glue("{M$chosen_path}/Analysis/cog_assignment_evidence.tsv"))
venom_evidence <- cog_evidence |> filter(Group == "venom_component")
venom_gos <- inner_join(venom_evidence, read_tsv(M$go_reference), by = join_by(x$Evidence == y$GO_IDs))

# ----------------------------------------
# Merging results
nd_merged_path <- glue("{M$wd}/results/ND_MERGED")
nd_run <- get_run(M$prefixes[[4]], M$ndpath)
ndm_run <- get_run(M$prefixes[[4]], nd_merged_path)

all_tests <- tibble()
stats <- tibble()
joined <- lapply(c("first", "second"), \(x) {
  tmp_joined <- inner_join(nd_run[[x]], ndm_run[[x]], by = join_by(header), suffix = c(".nd", ".ndm")) |>
    mutate(pass = x)
  test <- with(tmp_joined, wilcox.test(pcoverage_align.nd,
    pcoverage_align.ndm,
    paired = TRUE, alternative = "less"
  )) |>
    htest2tb(data.name = "ND x ND merged", alternative = "ND less") |>
    mutate(pass = x)
  all_tests <<- bind_rows(all_tests, test)
  tmp_stats <- with(tmp_joined, tibble(
    equals = sum(pcoverage_align.nd == pcoverage_align.ndm),
    nd_greater = sum(pcoverage_align.nd > pcoverage_align.ndm),
    ndm_greater = sum(pcoverage_align.nd < pcoverage_align.ndm),
    pass = x
  ))
  stats <<- bind_rows(stats, tmp_stats)

  tmp_joined
}) |>
  bind_rows()

stats <- mutate(stats, prop_ndm_greater = ndm_greater / (ndm_greater + equals + nd_greater)) |> gt()

TABLES$nd_ndm_comparison_stats <- stats

default_run <- get_run("C_indra", M$path)
default <- bind_rows(
  mutate(default_run$first, pass = "first"),
  mutate(default_run$second, pass = "second")
)


joined_def <- joined |>
  select(header, pcoverage_align.ndm) |>
  inner_join(default, by = join_by(header))

ndm_default_comparison <- compare_vals_x_y("default", "ND merged", "pcoverage_align",
  "pcoverage_align.ndm", joined_def,
  color = "pass",
  continuous = FALSE, segment_color = "black"
) + M$default_theme +
  theme(
    axis.title.y = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  ) +
  xlab("Coverage (%), default")

nd_ndm_comparison <- compare_vals_x_y(
  "ND", "ND merged", "pcoverage_align.nd",
  "pcoverage_align.ndm", joined,
  color = "pass", continuous = FALSE, segment_color = "black"
) +
  M$default_theme +
  xlab("Coverage (%), ND") +
  ylab("Coverage (%), ND merged") +
  guides(color = "none")

diff <- joined$pcoverage_align.ndm - joined$pcoverage_align.nd
summary(diff[diff != 0])

GRAPHS$nd_ndm_comparison <- cowplot::plot_grid(nd_ndm_comparison, ndm_default_comparison)
attr(GRAPHS$nd_ndm_comparison, "width") <- 17

# Show that peptides are lost in default
get_combined <- function(path) {
  lst <- list(
    first = read_tsv(glue("{path}/{M$passes[[1]]}/Combined/database_hits.tsv")),
    second = read_tsv(glue("{path}/{M$passes[[2]]}/Combined/database_hits.tsv"))
  )
  merge_runs(lst)
}

nd_all <- get_combined(M$ndpath) |> mutate(num_peps = str_count(peptideIds, ";") + 1)
def_all <- get_combined(M$path) |> mutate(num_peps = str_count(peptideIds, ";") + 1)


n_peps_compare <- inner_join(nd_all, def_all, by = join_by(header))

test <- with(n_peps_compare, wilcox.test(num_peps.x, num_peps.y,
  paired = TRUE,
  alternative = "greater"
)) |> htest2tb(data.name = "Peptide number ND x default", alternative = "ND greater")
TABLES$nd_d_peptide_number_test <- test

n_peps <- list(default = n_peps_compare$num_peps.y, ND = n_peps_compare$num_peps.x)

# ----------------------------------------
# Why is nd_merged so much better than default and has higher coverage than ND?
# 1. Could it be that the proteins that benefit from the increased coverage in
# nd_merged simply aren't present in default?
# 2. Does nd_merged improve the coverage of proteins that
# had higher coverage in ND than in default or is it other proteins?

better_ndm <- joined |>
  filter(pcoverage_align.ndm > pcoverage_align.nd) |>
  distinct(header, .keep_all = TRUE)
jd2 <- inner_join(nd_run$second, default_run$second, by = join_by(header)) |> distinct(ProteinId.x, .keep_all = TRUE)

nd_better <- jd2 |> filter(pcoverage_align.x > pcoverage_align.y)
d_better <- jd2 |> filter(pcoverage_align.y > pcoverage_align.x)

nd_better |>
  select(contains("Matched")) |>
  pluck("MatchedPeptideIds.x")

nd_better
# Answer to 2
improved_by_ndm <- nd_better |> filter(header %in% better_ndm$header)
already_better_in_nd <- d_better |> filter(header %in% better_ndm$header)
not <- nd_better |> filter(!header %in% better_ndm$header)

# Check open search peptides
# Make this into a function and avg across params
count_open_peps <- function(path, param) {
  peptide_map_all_path <- glue("{path}/{M$chosen_pass}/percolator_peptide_map_all.tsv")
  pmap_all <- read_tsv(peptide_map_all_path)
  open <- filter(pmap_all, engine %in% M$open_search_engines)
  standard <- filter(pmap_all, !engine %in% M$open_search_engines)
  n_open <- nrow(open)
  n_standard <- nrow(standard)
  new_peps <- open |> filter(!peptideIds %in% standard$peptideIds)
  n_all_prot <- open$ProteinId |>
    unique() |>
    length()
  new_prot <- open |>
    filter(!ProteinId %in% standard$ProteinId) |>
    pluck("ProteinId") |>
    unique()
  new_mod_peps <- open |> filter(!modifiedPeptideIds %in% standard$modifiedPeptideIds)
  n_new_peps <- nrow(new_peps)
  n_new_mod_peps <- nrow(new_mod_peps)
  tibble(
    `Run` = param, `N new peptides (unmodified)` = n_new_peps, `Proportion (unmodified)` = n_new_peps / n_open,
    `N new peptides` = n_new_mod_peps, `Proportion` = n_new_mod_peps / n_open,
    `Proportion of new proteins` = length(new_prot) / n_all_prot
  )
}

open_stats <- lapply(seq_len(length(M$all_paths)), \(x) {
  count_open_peps(M$all_paths[[x]], M$params[[x]])
}) |>
  bind_rows() |>
  mutate(across(where(is.double), \(x) round(x, 2)))

TABLES$open_stats <- open_stats |> gt()

# Check if de novo peptides and the proteins they match to go into the same protein groups

path <- M$chosen_path
tb <- read_tsv(glue("{path}/{M$chosen_pass}/percolator_all.tsv"))
mod_counts <- tb |>
  filter(!is.na(mods)) |>
  pluck("mods") |>
  lapply(\(x) str_split_1(x, ";")) |>
  unlist() |>
  map_chr(\(x) str_remove(x, "\\|[0-9]+$")) |>
  table()

mod_counts


save(TABLES, M$outdir)
