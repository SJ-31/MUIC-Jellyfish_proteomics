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


membrane <- read_lines(glue("{DATA_DIR}/membrane")) %>% paste0(., " membrane")
process <- read_lines(glue("{DATA_DIR}/processes")) %>% paste0(., " process")
system <- read_lines(glue("{DATA_DIR}/systems")) %>% paste0(., " system")
response <- read_lines(glue("{DATA_DIR}/response")) %>% paste0("response ", .)
pathway <- read_lines(glue("{DATA_DIR}/pathway")) %>% paste0(., " pathway")
cell <- c(
  read_lines(glue("{DATA_DIR}/cell_post")) %>% paste0("cell ", .),
  read_lines(glue("{DATA_DIR}/cell_pre")) %>% paste0(., " cell")
)
other <- c(
  "amino acid", "fatty acid",
  "cuticle chitin", "omega N", "gaseous exchange",
  "smooth muscle", "cardiac muscle",
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
  "secretory granule", "light stimulus",
  "external stimulus", "abiotic stimulus", "mechanical stimulus",
  "lytic vacuole", "pH reduction", "pigment granule",
  "intracellular vesicle", "pole plasm", "developmental growth"
)
PHRASES <- phrase(c(other, membrane, cell, process, system, response, pathway))
REPLACEMENTS <- quanteda::dictionary(
  list(
    cytotoxicity = "cellular cytotoxicity",
    senescence = "cellular senescence",
    cell_morphogenesis = c("cellular morphogenesis", "cell morphogenesis")
  )
)
UNWANTED <- c("cellular", "activity", "regulation", "of")

tbCompoundWord <- function(tb, word, fun) {
  has <- tb %>% filter(grepl(word, term))
  has_not <- tb %>% filter(!grepl(word, term))
  has$term <- fun(has$term)
  bind_rows(has, has_not)
}

modifyStr <- function(str) {
  assertArg(str, \(x) is.atomic(x))
  str %>%
    quanteda::corpus() %>%
    quanteda::tokens() %>%
    quanteda::tokens_compound(., pattern = PHRASES) %>%
    quanteda::tokens_lookup(., dictionary = REPLACEMENTS, exclusive = FALSE) %>%
    quanteda::tokens_remove(., pattern = PHRASES) %>%
    as.list() %>%
    map_chr(., \(x) str_c(x, collapse = " ")) %>%
    `names<-`(NULL)
}
