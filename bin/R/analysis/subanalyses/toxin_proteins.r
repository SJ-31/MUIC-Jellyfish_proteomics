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


rename_toxins <- swap_w_map(c(
  "c_type_lectin" = "C-type lectin", "nerve_growth_factor" = "nerve growth factor",
  "pesticidal_crystal" = "pesticidal crystal", "pore_forming" = "pore-forming",
  "serine_protease" = "serine protease", "amino_acid_oxidase" = "amino acid oxidase"
))
toxin_groups <- read_tsv(glue("{M$outdir}/toxin_groups.tsv")) |>
  filter(!is.na(Group)) |>
  mutate(Group = rename_toxins(Group), Group = map_chr(Group, str_to_title))


rename_tax <- swap_w_map(c("Cnidaria;Chytridiomycota" = "Cnidaria"))
w_intensity <- M$data |>
  inner_join(M$taxa_tb, by = join_by(ProteinId)) |>
  inner_join(merge_lfq(M$data, "mean"), by = join_by(ProteinId)) |>
  select(ProteinId, log_intensity, organism, header, GroupUP, all_of(M$taxa_cols)) |>
  filter(!is.na(log_intensity)) |>
  group_by(GroupUP) |>
  summarize(
    mean_intensity = mean(log_intensity, na.rm = TRUE), # Get mean
    # intensity over the entries in the group (since we don't know which single entry is
    # actually present
    organism = list(organism), header = list(header),
    across(all_of(M$taxa_cols), list)
  ) |>
  mutate(across(
    all_of(M$taxa_cols),
    \(x) lapply(x, \(y) rename_tax(modes_concat((y)))) |> unlist()
  ))

toxin_tb <- w_intensity |>
  inner_join(toxin_groups, by = join_by(GroupUP)) |>
  mutate(Group = map2_chr(Group, organism, \(g, o) {
    if (any(str_detect(o, "Chironex")) && g == "Pore Forming") {
      "Jft"
    } else {
      g
    }
  }))


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

tax_col <- "Phylum"
toxin_w_tax <- toxin_tb |>
  select(Group, {{ tax_col }}) |>
  group_by(Group) |>
  count(!!as.symbol(tax_col), name = "taxa_counts") |>
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
  arrange(desc(n)) |>
  pluck("Group") |>
  unique()

toxin_stacked_bar <- toxin_w_tax |>
  ggplot(aes(fill = !!as.symbol(tax_col), y = log_sum_intensity, x = Group)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits = group_order) +
  scale_fill_paletteer_d(PALETTE) +
  ylab("log intensity") +
  xlab("Toxin group")

GRAPHS$toxin_stacked_bar <- toxin_stacked_bar
attr(GRAPHS$toxin_stacked_bar, "width") <- 15

GRAPHS$toxin_tm_plotly <- toxin_tm

hemolysins <- toxin_tb |> filter(Group == "Hemolysin")
pore_forming <- toxin_tb |> filter(Group == "Pore Forming")

to_chord <- toxin_w_tax |>
  rename(from = Group, to = !!as.symbol(tax_col)) |>
  select(from, to)
# https://yjunechoe.github.io/posts/2020-06-30-treemap-with-ggplot/
# save_chord(to_chord,
#   filename = glue("{M$outdir}/taxonomy/toxin_chord.svg"), width = 15, height = 15
# )

py_plots$plotly_sunburst(
  to_sb, list(
    file = glue("{M$outdir}/taxonomy/toxin_sb_plotly.html"),
    width = 1500, height = 1200
  ),
  list(title = list(text = "Toxin intensities"))
)
py_plots$plotly_treemap(to_sb,
  list(
    file = glue("{M$outdir}/taxonomy/toxin_tm_plotly.html"),
    width = 1500, height = 1200
  ),
  layout_params = list(title = list(text = "Toxin intensities"))
)
save(c(TABLES, GRAPHS), glue("{M$outdir}/taxonomy"))
