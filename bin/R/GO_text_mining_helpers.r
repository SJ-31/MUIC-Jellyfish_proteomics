library(quanteda)
library(ggwordcloud)
library(tidytext)
DATA_DIR <- args$go_tm_dir

compoundSpecial <- function(str, pattern, expected_matches) {
  assertArg(str, \(x) is.atomic(x))
  match <- str_match(str, pattern)
  selected <- lapply(seq(2, expected_matches + 1), \(x)  {
    match[, x] %>%
      str_trim() %>%
      str_replace_all(., " ", "_")
  })
  purrr::reduce(selected, \(x, y) paste0(x, " ", y)) %>% str_trim()
}

sandwichText <- function(str, leftFun, middleFun, rightFun) {
  assertArg(str, \(x) is.atomic(x))
  extracted <- lapply(
    c(leftFun, middleFun, rightFun),
    \(F) F(str) %>% map_chr(., \(x) ifelse(is.na(x), "", x))
  )
  purrr::reduce(extracted, \(x, y) paste0(x, " ", y))
}

getText <- function(file_name) {
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
  lapply(getText) %>%
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
REPLACEMENTS <- quanteda::dictionary(
  list(
    cytotoxicity = "cellular cytotoxicity",
    senescence = "cellular senescence",
    cell_morphogenesis = c("cellular morphogenesis", "cell morphogenesis")
  )
)
UNWANTED <- c(
  "cellular", "activity", "regulation", "of", "process",
  "molecule", "substance", "primary", "nucleic", "non",
  "positive", "negative", "biosynthetic", "catabolic", "metabolic",
  "cell", "molecular", "active", "entity", "male", "multi",
  "membrane", "containing", "bound", "intracellular", ""
)

tbCompoundWord <- function(tb, word, fun) {
  has <- tb %>% filter(grepl(word, term))
  has_not <- tb %>% filter(!grepl(word, term))
  has$term <- fun(has$term)
  bind_rows(has, has_not)
}

compoundAll <- function(str) {
  assertArg(str, \(x) is.atomic(x))
  str %>%
    quanteda::corpus() %>%
    quanteda::tokens() %>%
    quanteda::tokens_compound(., pattern = PHRASES) %>%
    quanteda::tokens_lookup(., dictionary = REPLACEMENTS, exclusive = FALSE) %>%
    quanteda::tokens_remove(., pattern = UNWANTED) %>%
    as.list() %>%
    map_chr(., \(x) str_c(x, collapse = " ")) %>%
    `names<-`(NULL) %>%
    map_chr(., \(x) { # Get only words that have been compounded
      if (str_detect(x, "_")) {
        extracted <- str_extract(x, "(\\b[a-zA-Z_]*_[a-zA-Z_]*\\b)")
        return(extracted)
      } else {
        return(x)
      }
    })
}


# TODO Complete the abbreviation functions
# Only abbreviations present in the data should be added
ABBREVS <- c(
  "localization" = "loc.",
  "organization" = "org.",
  "establishment" = "est.",
  "processing" = "prcs.",
  "biological" = "biol.",
  "intracellular" = "intr.",
  "activity" = "act.",
  "extracellular" = "extr.",
  "developmental" = "dev.",
  "proliferation" = "prol.",
  "response" = "resp.",
  "migration" = "mig.",
  "maintenance" = "main.",
  "signaling" = "sign.",
  "regulation" = "reg.",
  "homeostasis" = "hom.",
  "differentiation" = "diff."
)
F_ABBREVS <- local({
  lst <- ABBREVS
  names(lst) <- names(ABBREVS) %>% map_chr(., \(x) glue("\\b{x}\\b"))
  lst
})

tokenize2Plot <- function(tb, params = NULL) {
  min_count <- lget(params, "min_count", quantile(tb$n, 0.80, na.rm = TRUE))
  top_n <- lget(params, "top_n", 50)
  abbreviate <- lget(params, "abbreviate", TRUE)
  sort_by <- lget(params, "sort_by", "n")
  term_col <- lget(params, "term_col", "term")
  compound <- lget(params, "compound", TRUE)
  result <- list()
  if (compound) {
    tb <- tb %>%
      mutate(., term = compoundAll(!!as.symbol(term_col))) %>%
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
    found_abbrevs <- tb$token %>%
      lapply(., \(x) str_split_1(x, " ")) %>%
      unlist() %>%
      keep(., \(x) x %in% names(ABBREVS))
    if (length(found_abbrevs) != 0) {
      legend_text <- found_abbrevs %>% map_chr(., \(x)  {
        key <- ABBREVS[x]
        glue("{key} = {x}")
      })
      result$abbrevs <- grid::legendGrob(legend_text)
    }
    tb$token <- tb$token %>% map_chr(., \(x) str_replace_all(x, F_ABBREVS))
  }
  result$tb <- tb
  return(result)
}

wordcloudCustom <- function(tb, params, abbrev_legend = NULL) {
  color_col <- lget(params, "color_col", "type")
  max_size <- lget(params, "max_size", 40)
  grid_size <- lget(params, "grid_size", 15)
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
    plot <- plot + guides(custom = guide_custom(abbrev_legend, title = "Abbreviations"))
  }
  plot
}


#' Higher-level interface for getting GO word clouds from results tbs
specialGoClouds <- function(
    tb, go_column = "GO_IDs",
    plot_params = list(),
    is_info_tb = FALSE,
    plot_special = TRUE, frequencies = NULL) {
  assertArg(tb, \(x) ("tbl" %in% class(x)) && (go_column %in% colnames(x)))
  gos <- goVector(tb, go_column = "GO_IDs")
  if (is.null(frequencies)) {
    frequencies <- gos %>%
      table() %>%
      table2Tb(., "GO_IDs")
  }
  if (!is_info_tb) {
    info_tb <- goInfoTb(unique(gos)) %>% inner_join(., frequencies)
  } else {
    info_tb <- inner_join(tb, frequencies)
  }
  to_plot <- list()
  # Plot special GO terms associated with a central concept, are
  # consistently described with adjectives e.g. "positive" or
  # "catabolic"
  bp <- info_tb %>% filter(ontology == "BP")

  to_plot$bp_regulation <- bp %>%
    filter(grepl("regulation", term)) %>%
    mutate(type = map_chr(term, \(x) {
      case_when(
        str_detect(x, "positive") ~ "positive",
        str_detect(x, "negative") ~ "negative",
        .default = "unspecified effect"
      )
    }))

  to_plot$bp_process <- bp %>%
    filter(grepl("process", term) & grepl("biosynthetic|metabolic|catabolic", term)) %>%
    mutate(type = map_chr(term, \(x) {
      case_when(
        str_detect(x, "catabolic") ~ "catabolic",
        str_detect(x, "biosynthetic") ~ "anabolic",
        .default = "unspecified metabolism"
      )
    }))
  tokens <- lapply(to_plot, \(x) tokenize2Plot(x, plot_params))
  to_plot$main <- local({
    all <- lapply(to_plot, \(x) x$GO_IDs) %>% unlist()
    info_tb %>% filter(!GO_IDs %in% all)
  })
  tokens$main <- tokenize2Plot(to_plot$main, list(compound = FALSE))
  clouds <- list()
  plot_params_copy <- plot_params
  plot_params_copy$color_col <- "type"
  clouds$bp_regulation <- wordcloudCustom(
    tokens$bp_regulation$tb,
    plot_params_copy,
    abbrev_legend = tokens$bp_regulation$abbrevs
  )
  clouds$bp_process <- wordcloudCustom(
    tokens$bp_process$tb,
    plot_params_copy,
    abbrev_legend = tokens$bp_process$abbrevs
  )
  plot_params$color_col <- "ontology"
  clouds$main <- wordcloudCustom(tokens$main$tb, plot_params,
    abbrev_legend = tokens$main$abbrevs
  )
  return(clouds)
}
