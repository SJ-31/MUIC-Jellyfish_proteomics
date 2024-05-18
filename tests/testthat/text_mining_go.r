library(quanteda.textplots)
library(quanteda.textstats)
library(quanteda)
library(tidytext)
library(simplifyEnrichment)
th <- new.env()
source(glue("{args$r_source}/GO_text_mining_helpers.r"), local = th)

frequencies <- goVector(d$sample_tb, go_column = "GO_IDs") %>%
  table() %>%
  table2Tb(., "GO_IDs")

sample <- d$go_vec$sample[1:1000]

info_tb <- goInfoTb(sample) %>% inner_join(., frequencies)



#'  Compound biologically meaningful phrases
compounded_tb <- info_tb %>% mutate(term = th$modifyStr(term))


special <- list()
# Special GO terms associated with a central concept, and may also be
# consistently described with adjectives e.g. "positive regulation of" or
# "catabolic"
bp <- info_tb %>% filter(ontology == "BP")
special$bp_regulation <- bp %>%
  filter(grepl("regulation", term)) %>%
  mutate(type = map_chr(term, \(x) {
    case_when(
      str_detect(x, "positive") ~ "positive",
      str_detect(x, "negative") ~ "negative",
      .default = "unspecified"
    )
  }))
special$bp_process <- bp %>%
  filter(grepl("process", term) & grepl("biosynthetic|metabolic|catabolic", term)) %>%
  mutate(type = map_chr(term, \(x) {
    case_when(
      str_detect(x, "catabolic") ~ "catabolic",
      str_detect(x, "biosynthetic") ~ "anabolic",
      .default = "unspecified"
    )
  }))

others <- local({
  all <- lapply(special, \(x) x$GO_IDs) %>% unlist()
  info_tb %>% filter(!GO_IDs %in% all)
})

# Need to filter those that are also related to "activity"
# special$bp_regulation <- local({
#   has_activity <- special$bp_regulation %>% filter(., grepl("activity", term))
#   others <- special$bp_regulation %>% filter(!GO_IDs %in% has_activity$GO_IDs)
#   has_activity$term <- th$compoundSpecial(
#     has_activity$term,
#     "^(?:.*of)?(.*)activity(?:.*in)?(.*)?", 2
#   )
#   others$term <- compoundSpecial(others$term, "^(?:.*of)?(.*)(?:.*from)?(.*)", 2)
#   bind_rows(others, has_activity)
# })

# special$bp_process$term <- th$compoundSpecial(
#   special$bp_process$term,
#   "^(?:.*of)?(.*)(?:metabolic|catabolic|biosynthetic.*)", 1
# )

bp_process_corpus <- corpus(special$bp_process,
  docid_field = "GO_IDs", text_field = "term"
)
docvars(bp_process_corpus) <- dplyr::select(special$bp_process, type)
bp_process_corpus %>%
  tokens() %>%
  dfm() %>%
  dfm_group(groups = type) %>%
  textplot_wordcloud(comparison = TRUE, )


corpus <- quanteda::corpus(others,
  docid_field = "GO_IDs",
  text_field = "term"
)
docvars(corpus) <- dplyr::select(info_tb, ontology, n)
tokens <- tokens(corpus) %>% tokens_compound(., pattern = PHRASES)
dfm <- dfm(tokens)

corpus %>%
  tokens(remove_punct = TRUE) %>%
  dfm() %>%
  dfm_group(groups = ontology) %>%
  textplot_wordcloud(comparison = TRUE)

# You'll need to use specific regexes and patterns to tokenize the definition
# or you could just use the term name?

# Create a word cloud based on the term names, adding specific key words
# onto the axes like "Regulation of" or "Response to"
# Something like this https://stackoverflow.com/questions/52804450/comparison-of-two-groups-using-word-cloud-comparison-r

# 2024-05-17 Plan:
# Word clouds for large go terms
# One word cloud for each subontology, which will display the most frequent
# nouns/phrases
# For BP, there will be two word clouds: one for regulation and response
# These will be comparative word clouds, positive - negative
# And if possible, catabolic vs biosynthetic (anabolic)

info_tb_tokenized <- info_tb %>% tidytext::unnest_tokens(., word, term)
info_tb_tokenized %>%
  anti_join(stop_words) %>%
  count(word, sort = FALSE)
