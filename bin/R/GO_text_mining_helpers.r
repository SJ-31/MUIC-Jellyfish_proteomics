library("quanteda")
library("ggplot2")
library("ggwordcloud")
library("tidytext")
DATA_DIR <- args$go_tm_dir

get_text <- function(file_name) {
  text <- read_lines(file_name)
  split <- str_split_1(text[1], "-")
  if (split[1] == "PRE") {
    return(paste0(text[-1], " ", split[2]))
  } else if (split[1] == "POST") {
    return(paste0(split[2], " ", text[-1]))
  } else {
    stop("Must specify reading mode on first line!")
  }
}

texts <- list.files(DATA_DIR, full.names = TRUE) %>%
  lapply(get_text) %>%
  unlist()

other <- c(
  "amino acid", "fatty acid",
  "cuticle chitin", "omega N", "gaseous exchange",
  "smooth muscle", "cardiac muscle",
  "plasma membrane bounded",
  "transmembrane transporter", "endoplasmic reticulum",
  "extracellular matrix", "spindle pole", "antigen processing",
  "DNA replication", "synaptic vesicle", "calcium ion",
  "transport vesicle", "viral capsid", "basal body",
  "myelin sheath",
  "mitotic spindle", "RNA polymerase", "DNA polymerase",
  "transcription factor", "vesicle lumen", "ubiquitin ligase",
  "actin filament", "receptor signaling",
  "growth factor", "skeletal muscle", "neural retina", "complex assembly",
  "cellular fusion", "Golgi apparatus", "visible light",
  "action potential", "synaptic signaling", "Golgi vesicle",
  "secretory granule", "light stimulus", "heart rate",
  "external stimulus", "abiotic stimulus", "mechanical stimulus",
  "lytic vacuole", "pH reduction", "pigment granule",
  "intracellular vesicle", "pole plasm", "developmental growth"
)
PHRASES <- quanteda::phrase(c(other, texts))

tb_compound_word <- function(tb, word, fun) {
  has <- tb %>% filter(grepl(word, term))
  has_not <- tb %>% filter(!grepl(word, term))
  has$term <- fun(has$term)
  bind_rows(has, has_not)
}

compound_all <- function(str) {
  check_arg(str, \(x) is.atomic(x))
  str_replace_all(str, " ", "_")
}


# TODO Complete the abbreviation functions
# Only abbreviations present in the data should be added
ABBREVS <- unlist(jsonlite::read(M$go_abbrev_file))
F_ABBREVS <- local({
  lst <- ABBREVS
  names(lst) <- names(ABBREVS) %>% map_chr(., \(x) glue("\\b{x}\\b"))
  lst
})

find_abbrev <- function(str) {
  lapply(names(ABBREVS), \(x) str_extract(str, x)) |>
    unlist() |>
    discard(is.na)
}

tokenize2plot <- function(tb, params = NULL) {
  min_count <- lget(params, "min_count", quantile(tb$n, 0.80, na.rm = TRUE))
  top_n <- lget(params, "top_n", 50)
  abbreviate <- lget(params, "abbreviate", TRUE)
  sort_by <- lget(params, "sort_by", "n")
  term_col <- lget(params, "term_col", "term")
  compound <- lget(params, "compound", TRUE)
  result <- list()
  if (compound) {
    tb <- tb %>%
      mutate(., term = compound_all(!!as.symbol(term_col))) %>%
      tidytext::unnest_tokens(., token, term) %>%
      filter(., !token %in% UNWANTED & nchar(token) > 2) %>%
      mutate(token = map_chr(token, \(x) str_replace_all(x, "_", " ")))
  } else {
    tb <- dplyr::rename(tb, token = !!as.symbol(term_col)) %>% filter(., !is.na(token))
  }
  tb <- tb %>%
    {
      if (sort_by == "n") filter(., n > min_count) else .
    } %>%
    distinct(token, .keep_all = TRUE) %>%
    arrange(., desc(!!as.symbol(sort_by))) %>%
    slice(1:top_n)
  if (abbreviate) {
    found_abbrevs <- tb$token |>
      lapply(find_abbrev) |>
      unlist() |>
      discard(is.na)
    if (length(found_abbrevs) != 0) {
      legend_text <- found_abbrevs %>%
        map_chr(., \(x)  {
          key <- ABBREVS[x]
          glue("{key} = {x}")
        }) %>%
        unique()
      result$abbrevs <- grid::legendGrob(legend_text)
    }
    tb$token <- tb$token %>% map_chr(., \(x) str_replace_all(x, F_ABBREVS))
  }
  result$tb <- tb
  return(result)
}

wordcloud_custom <- function(tb, params, abbrev_legend = NULL) {
  color_col <- lget(params, "color_col", "type")
  max_size <- lget(params, "max_size", 40)
  grid_size <- lget(params, "grid_size", 25)
  max_grid_size <- lget(params, "max_grid_size", 200)
  size <- lget(params, "word_size", "n")
  shape <- lget(params, "shape", "pentagon")
  if (color_col %in% colnames(tb)) {
    plot <- ggplot(tb, aes(
      angle_group = !!as.symbol(color_col),
      label = token, size = !!as.symbol(size),
      color = !!as.symbol(color_col)
    ))
  } else {
    plot <- ggplot(tb, aes(
      label = token, size = !!as.symbol(size)
    ))
  }
  plot <- plot +
    ggwordcloud::geom_text_wordcloud_area(
      show.legend = TRUE, family = "serif",
      grid_size = grid_size,
      max_grid_size = max_grid_size,
      rm_outside = TRUE,
      shape = shape
    ) +
    theme(panel.background = element_rect(fill = "#eff1f5")) +
    scale_size_area(max_size = max_size) +
    guides(size = FALSE, alpha = FALSE) +
    theme(text = element_text(family = "Ubuntu"))
  if (!is.null(abbrev_legend)) {
    plot <- plot + guides(custom = ggplot2::guide_custom(abbrev_legend, title = "Abbreviations"))
  }
  plot
}


LEVEL_MAP <- glue("{M$wd}/data/reference/go_level_map.tsv")
MIN_MAX_LEVEL <- c(3, 7) # Closed interval for GO levels
MAP_REF_FILE <- glue("{M$wd}/data/reference/go_map_generic.tsv")
#' Higher-level interface for getting GO word clouds from results tbs
special_go_clouds <- function(
    tb, go_column = "GO_IDs",
    plot_params = list(),
    is_info_tb = FALSE,
    plot_special = TRUE, frequencies = NULL, to_generic = TRUE) {
  check_arg(tb, \(x) ("tbl" %in% class(x)) && (go_column %in% colnames(x)))
  gos <- get_go_vec(tb, go_column = "GO_IDs") |> discard(\(x) x %in% UNWANTED)
  if (to_generic) {
    ref <- read_tsv(MAP_REF_FILE)
    map <- hash::hash(key = ref["child"], values = ref["parent"])
    gos <- map_chr(gos, \(x) ifelse(hash::has.key(x, map), map[[x]], x))
  }
  if (is.null(frequencies)) {
    frequencies <- gos %>%
      table() %>%
      table2tb(., "GO_IDs")
  }
  if (exists("args") && !is.null(args$go_info)) {
    info_tb <- read_tsv(args$go_info) |>
      filter(GO_IDs %in% gos) |>
      inner_join(frequencies)
  } else if (!is_info_tb) {
    info_tb <- go_info_tb(unique(gos)) %>% inner_join(., frequencies)
  } else {
    info_tb <- inner_join(tb, frequencies)
  }
  info_tb <- inner_join(info_tb, read_tsv(LEVEL_MAP), by = join_by(GO_IDs)) |>
    filter((level >= MIN_MAX_LEVEL[1]) & (level <= MIN_MAX_LEVEL[2]))
  gos <- gos |> keep(\(x) x %in% info_tb$GO_IDs)
  to_plot <- list()
  # Plot special GO terms associated with a central concept, are
  # consistently described with adjectives e.g. "positive" or
  # "catabolic"
  bp <- info_tb %>% filter(ontology == "BP")

  to_plot$bp_regulation <- bp %>%
    filter(grepl("regulation", term)) %>%
    mutate(
      type = map_chr(term, \(x) {
        case_when(
          str_detect(x, "positive") ~ "positive",
          str_detect(x, "negative") ~ "negative",
          .default = "unspecified effect"
        )
      })
    )

  to_plot$bp_process <- bp %>%
    filter(grepl("process", term) & grepl("biosynthetic|metabolic|catabolic", term)) %>%
    mutate(type = map_chr(term, \(x) {
      case_when(
        str_detect(x, "catabolic") ~ "catabolic",
        str_detect(x, "biosynthetic") ~ "anabolic",
        .default = "unspecified metabolism"
      )
    }))
  tokens <- lapply(to_plot, \(x) tokenize2plot(x, plot_params))
  to_plot$main <- local({
    all <- lapply(to_plot, \(x) x$GO_IDs) %>% unlist()
    info_tb %>% filter(!GO_IDs %in% all)
  })
  tokens$main <- tokenize2plot(to_plot$main, list(compound = FALSE))
  clouds <- list()
  plot_params_copy <- plot_params
  plot_params_copy$color_col <- "type"
  clouds$bp_regulation <- wordcloud_custom(
    tokens$bp_regulation$tb,
    plot_params_copy,
    abbrev_legend = tokens$bp_regulation$abbrevs
  )
  clouds$bp_process <- wordcloud_custom(
    tokens$bp_process$tb,
    plot_params_copy,
    abbrev_legend = tokens$bp_process$abbrevs
  )
  plot_params$color_col <- "ontology"
  clouds$main <- wordcloud_custom(tokens$main$tb, plot_params,
    abbrev_legend = tokens$main$abbrevs
  )
  return(clouds)
}
