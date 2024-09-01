if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
py_plots <- new.env()
reticulate::source_python(glue("{M$python_source}/plotting.py"), envir = py_plots)

PALETTE <- "ggthemes::Classic_20"
TABLES <- list()
GRAPHS <- list()
# ----------------------------------------
# Toxin protein analysis

specify_other <- function(tb) {
  mutate(tb, Group = case_when(
    str_detect(str_to_lower(entry_name), "cftx|catx|cqtx|tx-like") ~ "JFT",
    str_detect(str_to_lower(entry_name), "ribonuclease") ~ "Ribonuclease",
    str_detect(str_to_lower(entry_name), "shkt") ~ "Metalloprotease",
    str_detect(str_to_lower(entry_name), "rtx toxin") ~ "Cytotoxin",
    str_detect(str_to_lower(entry_name), "peptidase c80") ~ "Cytotoxin",
    .default = Group
  ))
}


rename_toxins <- swap_w_map(c(
  "c_type_lectin" = "C-type lectin", "nerve_growth_factor" = "nerve growth factor",
  "pesticidal_crystal" = "pesticidal crystal", "pore_forming" = "pore-forming",
  "serine_protease" = "serine protease", "amino_acid_oxidase" = "amino acid oxidase"
))
toxin_groups <- read_tsv(glue("{M$chosen_path}/Analysis/toxin_groups.tsv")) |>
  filter(!is.na(Group)) |>
  mutate(Group = rename_toxins(Group), Group = map_chr(Group, str_to_title))

w_intensity |> filter(grepl("Trypsin", PFAMs)) |> see()

data <- read_tsv(M$data_w_cat_path) |>
  inner_join(M$lfq, by = join_by(ProteinId))

w_intensity <- data |>
  inner_join(M$taxa_tb, by = join_by(ProteinId)) |>
  inner_join(merge_lfq(data, "mean"), by = join_by(ProteinId)) |>
  select(ProteinId, log_intensity, organism, header, PFAMs, entry_name, GroupUP, all_of(M$taxa_cols)) |>
  filter(!is.na(log_intensity)) |>
  group_by(GroupUP) |>
  summarize(
    mean_intensity = mean(log_intensity, na.rm = TRUE), # Get mean
    # intensity over the entries in the group (since we don't know which single entry is
    # actually present
    organism = list(organism), header = list(header),
    entry_name = paste0(entry_name, collapse = ";"),
    PFAMs = paste0(PFAMs, collapse = ";"),
    across(all_of(M$taxa_cols), list)
  ) |>
  mutate(
    across(
      all_of(M$taxa_cols),
      \(x) lapply(x, \(y) modes(y, first = TRUE)) |> unlist()
    ),
    rank = rank(-mean_intensity),
    PFAMs = map_chr(PFAMs, split_unique_join)
  )


toxin_tb <- w_intensity |>
  inner_join(toxin_groups, by = join_by(GroupUP)) |>
  specify_other()

tox_w_other <- w_intensity |>
  full_join(toxin_groups, by = join_by(GroupUP)) |>
  specify_other()

tax_col <- "Phylum"

toxin_w_lfq <- toxin_tb |>
  group_by(Group, !!as.symbol(tax_col)) |>
  summarize(
    sum_intensity = sum(mean_intensity),
    mean_intensity = mean(mean_intensity),
  ) |>
  mutate(percent_intensity = sum_intensity / sum(sum_intensity)) |>
  group_by(Group) |>
  mutate(group_sum_intensity = sum(sum_intensity)) |>
  ungroup() |>
  mutate(log_sum_intensity = round(log(group_sum_intensity) * percent_intensity, 2))
# Since we group on both columns, the sum will based on `tax_col`
# adding up `tax_col` for a specific toxin group gives the `group_sum_intensity`
# of that toxin group

toxin_counts <- table(toxin_tb$Group) |>
  table2tb(id_col = "Group") |>
  inner_join(toxin_w_lfq)

toxin_w_tax <- toxin_tb |>
  select(Group, {{ tax_col }}) |>
  group_by(Group) |>
  dplyr::count(!!as.symbol(tax_col), name = "taxa_counts") |>
  ungroup() |>
  inner_join(toxin_counts) |>
  mutate(
    n_percent = n / sum(n),
    taxa_counts_percent = taxa_counts / sum(taxa_counts),
    #  log_sum_intensity
    # is the log of the summed intensity of entries belonging to a specific
    # `tax_col`, within a specific toxin group
  )

to_sb <- plotly_treemap_format(toxin_w_tax,
  outer_col = "Group",
  inner_col = tax_col, value_col = "sum_intensity"
)

group_order <- toxin_w_tax |>
  filter(Group != "Other") |>
  group_by(Group) |>
  summarise(log_sum_intensity = sum(log_sum_intensity)) |>
  arrange(desc(log_sum_intensity)) |>
  pluck("Group") |>
  unique()

toxin_stacked_bar <- toxin_w_tax |>
  filter(!is.na(Phylum)) |>
  ggplot(aes(fill = !!as.symbol(tax_col), y = log_sum_intensity, x = Group)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits = group_order) +
  scale_fill_paletteer_d(PALETTE) +
  ylab("log intensity") +
  theme(plot.margin = margin(1, 0, 0, 0, "cm")) +
  xlab("Toxin group")

toxin_rank_plot <- tox_w_other |>
  filter(!is.na(Group) & !is.na(!!as.symbol(tax_col))) |>
  ggplot(aes(x = Group, y = rank, color = !!as.symbol(tax_col))) +
  geom_point(size = 5, shape = 18) +
  scale_x_discrete(position = "top", limits = group_order) +
  scale_color_paletteer_d(PALETTE) +
  theme_bw() +
  M$default_theme +
  theme(
    axis.title.x = element_blank(), axis.text.x = element_blank(),
  ) +
  ylab("Rank (1 means highest intensity)") +
  guides(color = "none")

GRAPHS$all_tox <- cowplot::plot_grid(toxin_stacked_bar, toxin_rank_plot,
  ncol = 1, align = "v",
  labels = "AUTO", label_size = 15
)
attr(GRAPHS$all_tox, "width") <- 19
attr(GRAPHS$all_tox, "height") <- 12

GRAPHS$toxin_stacked_bar <- toxin_stacked_bar + M$default_theme

GRAPHS$toxin_stacked_bar
attr(GRAPHS$toxin_stacked_bar, "width") <- 19


sunburst <- py_plots$plotly_sunburst(
  to_sb, list(title = list(text = "Toxin intensities"))
)
set_attrs(sunburst, list(height = 1200, width = 1500))

GRAPHS$sunburst <- sunburst

GRAPHS$treemap <- py_plots$plotly_treemap(to_sb,
  layout_params = list(title = list(text = "Toxin intensities"))
)
set_attrs(GRAPHS$treemap, list(height = 1200, width = 1500))

save(c(TABLES, GRAPHS), glue("{M$outdir}/taxonomy"))

# EXTRA
grouped_cat <- read_tsv(glue("{M$chosen_path}/Analysis/grouped_annotations_only.tsv"))
go_data <- read_tsv(M$go_reference)
evidence <- read_tsv(glue("{M$chosen_path}/Analysis/toxin_groups_evidence.tsv"))
if (!"term" %in% colnames(evidence)) {
  evidence |>
    left_join(go_data, by = join_by(x$Evidence == y$GO_IDs)) |>
    write_tsv(glue("{M$chosen_path}/Analysis/toxin_groups_evidence.tsv"))
}

with_go <- function(tb) {
  separate_longer_delim(tb, GO_IDs, ";") |>
    group_by(GO_IDs) |>
    summarise(count = n()) |>
    inner_join(go_data, by = join_by(GO_IDs)) |>
    arrange(desc(count))
}

mito_groups <- grouped_cat |>
  filter(str_detect(str_to_lower(entry_name), "mitochondrial")) |>
  pluck("GroupUP")

joined_tox <- inner_join(toxin_groups, grouped_cat, by = join_by(GroupUP)) |>
  left_join(select(w_intensity, c(GroupUP, mean_intensity)), by = join_by(GroupUP)) |>
  mutate(intensity_percent = replace_na(mean_intensity / sum(mean_intensity, na.rm = TRUE), 0)) |>
  select(-c(interpro_accession, interpro_db, interpro_pathways)) |>
  specify_other() |>
  relocate(c(entry_name, pcoverage_align, intensity_percent, mean_intensity), .before = everything())

write_tsv(joined_tox, glue("{M$chosen_path}/Analysis/grouped_tox.tsv"))

get_sp_name <- function(str, genus_only = FALSE) {
  splits <- str_split_1(str, " ") |> discard(\(x) str_detect(x, "uncultured|sp\\.|strain"))
  if (genus_only) {
    return(splits[1])
  }
  if (length(splits) > 2) {
    glue("{splits[1]} {splits[2]}")
  } else {
    splits[1]
  }
}

# Check outer membrane proteins
ompa <- joined_tox |> filter(str_detect(entry_name, "Outer membrane") |
  str_detect(entry_name, "OmpA"))

ompa_genera <- flatten_by(ompa$organism, ";") |>
  map_chr(\(x) get_sp_name(x, genus_only = TRUE)) |>
  table()

ompa_species <- flatten_by(ompa$organism, ";") |>
  map_chr(get_sp_name) |>
  table()

w_cat <- read_tsv(M$data_w_cat_path)
lfq <- merge_lfq(inner_join(w_cat, M$lfq, by = join_by(ProteinId)), "mean")
tb <- w_cat %>%
  inner_join(lfq, by = join_by(ProteinId)) |>
  group_by(GroupUP) |>
  summarise(
    assigned_COG = modes(assigned_COG),
    header = paste0(header, collapse = ";"),
    size = n(),
    across(where(is.double), \(x) mean(x, na.rm = TRUE))
  ) |>
  filter(!is.na(log_intensity)) %>%
  arrange(log_intensity) %>%
  mutate(rank = seq_len(nrow(.)))

mito <- tb |> mutate(type = str_detect(header, "mitochondria"))
mito |> ggplot(aes(x = rank, y = log_intensity, color = type)) +
  geom_point()
