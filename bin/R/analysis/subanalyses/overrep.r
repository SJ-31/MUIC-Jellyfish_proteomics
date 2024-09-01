library("paletteer")
library("seqinr")
library("ggVennDiagram")
if (!exists("SOURCED")) {
  source(paste0(dirname(getwd()), "/", "all_analyses.r"))
  SOURCED <- TRUE
}
source(glue("{args$r_source}/GO_text_mining_helpers.r"))
source(glue("{args$r_source}/GO_chord.r"))
TABLES <- list()
GRAPHS <- list()

# ----------------------------------------
# Analyses for intensity-enriched terms
intensities <- c("high", "medium", "low")
intensity_tbs <- lapply(
  intensities,
  \(x) {
    read_tsv(glue("{M$ontologizer_path}/{x}_intensity.tsv")) |>
      mutate(subset = glue("{x} intensity"))
  }
) |>
  `names<-`(intensities)

intensity_vecs <- intensity_tbs |>
  lapply(\(x) {
    x |>
      filter(p.adjusted < 0.05) |>
      pluck("ID")
  })
data <- read_tsv(M$data_w_cat_path)

go_data <- read_tsv(M$go_reference)
ontologizer <- get_ontologizer(M$ontologizer_path)


mod_names <- c("Met_ox", "Nterm_acetyl", "Lys_acetyl")
enriched_mods <- lapply(mod_names, \(x) {
  read_tsv(glue("{M$ontologizer_path}/{x}.tsv")) |>
    mutate(subset = str_replace(x, "_", " "))
})

all_ontologizer <- bind_rows(
  mutate(ontologizer$id_with_open, subset = "id with open"),
  mutate(ontologizer$unknown_to_db, subset = "not DBP"),
  mutate(ontologizer$transcriptome, subset = "transcriptome"),
  mutate(ontologizer$denovo, subset = "denovo"),
) |>
  bind_rows(
    bind_rows(intensity_tbs),
    bind_rows(enriched_mods)
  ) |>
  inner_join(go_data, by = join_by(x$ID == y$GO_IDs)) |>
  filter(p.adjusted < 0.05) |>
  mutate(sorted_p = invert_p_values(p.adjusted))

intensity_vecs_ontology <- lapply(intensity_vecs, \(x) {
  ids_into_ontology(x, target = "GOID", collapse = FALSE)
})

ta <- new.env()
reticulate::source_python(glue("{args$python_source}/trace_alignments.py"), envir = ta)

subsets <- unique(all_ontologizer$subset)
id_lists <- lapply(subsets, \(x) filter(all_ontologizer, subset == x)$ID) |> `names<-`(subsets)


id2term <- with(go_data, setNames(term, GO_IDs))
slims <- get_go_slim(unique(all_ontologizer$ID), M$go_path, M$go_slim_path)
all_ontologizer <- mutate(all_ontologizer,
  slim = map_chr(all_ontologizer$ID, \(x) {
    find_slims <- slims[[x]]$all |> discard(\(y) y == x)
    if (length(find_slims) > 0) {
      base::sample(find_slims, size = 1)
    } else {
      "NONE"
    }
  }),
  slim_name = map_chr(slim, \(x) id2term[x]),
  ontology_full = case_match(
    ontology,
    "MF" ~ "Molecular Function",
    "CC" ~ "Cellular Component",
    "BP" ~ "Biological Process"
  )
)
# ----------------------------------------
# Check if each of the subsets are associated with one another in terms
# of protein groups

sort_pairs_overlap <- function(tb) {
  tb |>
    distinct(pair, .keep_all = TRUE) |>
    mutate(
      pair = mapply(\(x, y) paste0(str_sort(c(x, y)), collapse = "_"), first, second),
      first = gsub("_.*", "", pair),
      second = gsub(".*_", "", pair)
    )
}

id_overlap <- ta$calculate_overlaps(id_lists) |> as_tibble()
mod_names2 <- map_chr(mod_names, \(x) str_replace(x, "_", " "))
GRAPHS$subset_overlap <- id_overlap |>
  sort_pairs_overlap() |>
  filter(first %in% mod_names2 & second %in% mod_names2) |>
  ggplot(aes(x = first, y = second, fill = overlap)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_paletteer_c("ggthemes::Classic Area Red", name = "Overlap coefficient\n (Enriched terms)") +
  M$default_theme +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "inside",
    legend.position.inside = c(.70, .25),
    panel.grid = element_blank(), axis.title.x = element_blank()
  ) +
  scale_x_discrete(position = "top") +
  ylab("Subset")

GRAPHS$subset_overlap

mod_binary <- read_tsv(glue("{M$ontologizer_path}/mod_binary.tsv"))
colnames(mod_binary) <- colnames(mod_binary) |> map_chr(\(x) str_replace(x, "_", " "))

lfq <- merge_lfq(inner_join(data, M$lfq), "mean") |>
  inner_join(select(data, ProteinId, GroupUP), by = join_by(ProteinId)) |>
  group_by(GroupUP) |>
  summarize(
    log_intensity = mean(log_intensity, na.rm = TRUE),
  ) |>
  mutate(
    `high intensity` = log_intensity >= quantile(log_intensity, 0.75),
    `medium intensity` = log_intensity > quantile(log_intensity, 0.25) & log_intensity < quantile(log_intensity, 0.75),
    `low intensity` = ifelse(!`high intensity` & !`medium intensity`, TRUE, FALSE)
  )
paste_unique <- function(x) {
  paste0(unique(x), collapse = ";")
}

# Cell 1: Met ox, no Lys acetyl, Cell 2: Met ox, Lys acetyl
# Cell 3: No Met ox, no lys aceytl, Cell 4: No Met, Lys acetyl
ox_acetyl <- table(!mod_binary$`Met ox`, !mod_binary$`Lys acetyl`)
chi <- chisq.test(ox_acetyl) |> htest2tb(data.name = "Association between Met ox and Lys acetyl")
# Was significant

grouped <- data |>
  group_by(GroupUP) |>
  mutate(ProteinId = str_sub(ProteinId, 1, 1)) |>
  summarise(
    GO_IDs = paste0(GO_IDs, collapse = ";"),
    ID_method = paste_unique(ID_method),
    ProteinId = paste_unique(ProteinId)
  ) |>
  filter(!is.na(GO_IDs)) |>
  mutate(
    `not DBP` = ifelse(str_detect(ProteinId, "P"), FALSE, TRUE),
    `id with open` = ifelse(str_detect(ID_method, "standard"), FALSE, TRUE)
  )

subset_binary <- grouped |>
  left_join(lfq, by = join_by(GroupUP)) |>
  left_join(mod_binary, by = join_by(GroupUP)) |>
  select(-c(ID_method, GO_IDs, log_intensity, ProteinId)) |>
  replace_na_all(FALSE)

investigate_association <- subsets[subsets %in% colnames(subset_binary)]

combos <- combn(investigate_association, 2)

independence <- lapply(seq_len(ncol(combos)), \(i) {
  x <- combos[1, i]
  y <- combos[2, i]
  contigency <- table(subset_binary[[x]], subset_binary[[y]])
  if (length(contigency) != 4) {
    return(tibble())
  }
  chisq.test(contigency) |>
    htest2tb(data.name = glue("{x} x {y}")) |>
    mutate(left = x, right = y)
}) |>
  bind_rows() |>
  select(-c(alternative, null)) |>
  get_adjusted_p() |>
  mutate(method = "Chi-squared")

TABLES$subset_chi <- independence

odds <- independence |>
  filter(significant == 1) |>
  apply(1, \(x) {
    left <- x[["left"]]
    right <- x[["right"]]
    tab <- table(subset_binary[[left]], !subset_binary[[right]],
      dnn = c(left, right)
    )
    or <- get_odds_ratio(tab)
    tibble(
      pair = glue("{left} x {right}"), OR = or,
      interpretation = glue("Odds of '{right}' are {round(or, 2)} times higher in 'not {left}' than '{left}'")
    )
  }) |>
  bind_rows()

TABLES$subset_odds <- odds

ta <- new.env()
reticulate::source_python(glue("{M$python_source}/trace_alignments.py"), envir = ta)
group_sets <- lapply(subsets, \(s) {
  if (s %in% colnames(subset_binary)) {
    subset_binary |>
      filter(!!as.symbol(s)) |>
      pluck("GroupUP")
  } else {
    NULL
  }
}) |>
  `names<-`(subsets) |>
  discard(is.null)

group_overlap <- ta$calculate_overlaps(group_sets) |> as_tibble()
GRAPHS$group_overlap <- group_overlap |>
  sort_pairs_overlap() |>
  ggplot(aes(x = second, y = first, fill = overlap)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_paletteer_c("ggthemes::Red-Gold", name = "Overlap coefficient\n (Protein groups)") +
  M$default_theme +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "inside",
    legend.position.inside = c(.35, .70),
    panel.grid = element_blank(), axis.title.x = element_blank()
  ) +
  scale_y_discrete(position = "right") +
  ylab("Subset")


joined_overlap <- group_overlap |>
  sort_pairs_overlap() |>
  inner_join(sort_pairs_overlap(id_overlap), by = join_by(pair))

ggplot(joined_overlap, aes(x = overlap.x, y = overlap.y)) +
  geom_point() +
  xlab("Overlap coefficient (Group)") +
  ylab("Overlap coefficient (Enriched terms)")
# with(joined_overlap, wilcox.test(overlap.x, overlap.y, paired = TRUE))

# ----------------------------------------
# Plot word clouds

wc <- new.env()
reticulate::source_python(glue("{M$python_source}/word_clouds.py"), envir = wc)

format_for_wc <- function(tb, val_col = "sorted_p", max = 30) {
  tokenized <- tb |>
    group_by(name) |>
    summarise(
      sorted_p = mean(sorted_p), n = n(),
      ontology = dplyr::first(ontology),
      ontology_full = dplyr::first(ontology_full),
      subset = paste0(subset, collapse = ";")
    ) |>
    tokenize2plot(tokenize_params, term_col = "name", sort_by = val_col)
  tokenized$tb <- tokenized$tb |>
    arrange(desc(val_col)) |>
    slice(1:max)
  list(
    tokens = tb2named_list(tokenized$tb, "token", val_col),
    abbrevs = tb2named_list(tokenized$legend_text, "abbrev", "text"),
    tb = tokenized$tb,
    map = tokenized$tb |> tb2named_list("token", "ontology_full")
  )
}

# PTMs
onto_colors <- list(
  `Biological Process` = "dark:salmon_r",
  `Cellular Component` = "dark:#008000_r",
  `Molecular Function` = "dark:#0000FF_r"
)


get_wc <- function(tb, cmap_legend = FALSE) {
  prep <- tb |> format_for_wc()
  params <- list(
    category2colormap = onto_colors, item2category = prep$map,
    abbrev_size = 20, cmap_legend = cmap_legend
  )
  wc$word_cloud_main(prep$tokens, prep$abbrevs,
    params = params
  )
}

GRAPHS$denovo_wc <- get_wc(all_ontologizer |> filter(subset == "denovo"), TRUE)
GRAPHS$t_wc <- get_wc(all_ontologizer |> filter(subset == "transcriptome"))

# GRAPHS$denovo_wc <- get_wc()

# ptm_wc <- all_ontologizer |>
#   filter(subset %in% mod_names2) |>
#   format_for_wc()

# wc_params <- list(
#   category2colormap = onto_colors, item2category = ptm_wc$map,
#   abbrev_size = 20
# )
# GRAPHS$ptm_wc <- wc$word_cloud_main(
#   ptm_wc$tokens,
#   ptm_wc$abbrevs,
#   params = wc_params
# )

# # Identification method
# # Change this to be de novo and transcriptome peptides
# id_wc <- all_ontologizer |>
#   filter(subset %in% c("id with open", "not DBP")) |>
#   format_for_wc()
# wc_params3 <- list(
#   category2colormap = onto_colors, item2category = id_wc$map,
#   abbrev_size = 20, cmap_legend = FALSE
# )
# GRAPHS$id_wc <- wc$word_cloud_main(id_wc$tokens, id_wc$abbrevs,
#   params = wc_params3
# )

# # Intensity
# intensity_wc <- all_ontologizer |>
#   filter(grepl("intensity", subset)) |>
#   format_for_wc()

# wc_params2 <- list(
#   category2colormap = onto_colors, item2category = intensity_wc$map,
#   abbrev_size = 20, cmap_legend = FALSE
# )
# GRAPHS$intensity_wc <- wc$word_cloud_main(intensity_wc$tokens, intensity_wc$abbrevs, params = wc_params2)

# EXTRA

TABLES$all_ontologizer_sig <- all_ontologizer

shared_terms <- function(a, b) {
  a <- select(a, ID, term, subset)
  b <- select(b, ID, term, subset)
  full_join(a, b, by = join_by(ID))
}

transcriptome <- all_ontologizer |> filter(subset == "transcriptome")
dnp <- all_ontologizer |> filter(subset == "denovo")
shared <- shared_terms(transcriptome, dnp)

# Targets
get_targets <- list(
  D = list(
    serine_protease = c("GO:0017171", "GO:0004252", "GO:0006508"),
    intermediate_filament = c("GO:0045103", "GO:0045104", "GO:0005882"),
    collagenase = c("GO:0022617", "GO:0030574"),
    antioxidant = c(
      "GO:1990748", # cellular detoxification
      "GO:0016209" # antioxidant activity
    ),
    hemostasis = c("GO:1900047", "GO:0030195")
  ),
  T = list(
    actin = c("GO:0030029", "GO:0015629", "GO:0098862"),
    localization = c("GO:0060341", "GO:0051649"),
    transporter = c(
      "GO:0046907", # intracellular transport
      "GO:0098660" # inorganic ion transmembrane transport
    ),
    signal_transduction = c("GO:0035556", "GO:0005102"),
    cytochrome = c("GO:0070069")
  )
)


seqs <- read_tsv(M$data_w_cat_path) |>
  select(ProteinId, header, GO_IDs, seq) |>
  filter(grepl("D|T", ProteinId) & !is.na(GO_IDs) & !is.na(seq)) |>
  separate_longer_delim(GO_IDs, ";")

outdir <- glue("{M$outdir}/query_seqs")
for (type in names(get_targets)) {
  cur_spec <- get_targets[[type]]
  for (spec in names(cur_spec)) {
    filename <- glue("{outdir}/{type}_{spec}.fasta")
    go_vec <- cur_spec[[spec]]
    filtered <- seqs |>
      filter(grepl(type, ProteinId) & GO_IDs %in% go_vec) |>
      distinct(header, .keep_all = TRUE)
    seq_vec <- filtered$seq
    print(seq_vec)
    id_vec <- filtered$ProteinId
    write.fasta(as.list(seq_vec), id_vec, file.out = filename, as.string = TRUE)
  }
}



GRAPHS$intensity_overlap <- ggVennDiagram(intensity_vecs) + scale_fill_paletteer_c("ggthemes::Classic Red")

save(c(TABLES, GRAPHS), glue("{M$outdir}/misc"))
