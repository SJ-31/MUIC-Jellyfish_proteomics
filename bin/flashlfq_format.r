library(glue)
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dirname <- dirname(script.name)
source(glue("{script.dirname}/read_engines.r"))

parser <- OptionParser()
parser <- add_option(parser, c("--metamorpheus"),
  type = "character",
  help = "metamorpheus psms"
)
parser <- add_option(parser, c("--identipy"),
  type = "character",
  help = "identipy psms"
)
parser <- add_option(parser, c("--maxquant"),
  type = "character",
  help = "maxquant psms"
)
parser <- add_option(parser, c("--msfragger"),
  type = "character",
  help = "msfragger psms"
)
parser <- add_option(parser, c("--comet"),
  type = "character",
  help = "comet psms"
)
parser <- add_option(parser, c("--tide"),
  type = "character",
  help = "tide psms"
)
parser <- add_option(parser, c("-o", "--output"),
  type = "character",
  help = "Output file name"
)
parser <- add_option(parser, c("-m", "--msms_mapping"),
  type = "character",
  help = "MsMs mapping file"
)
args <- parse_args(parser)

mapping <- read.delim(args$msms_mapping, sep = "\t")


file_list <- list(comet = args$comet, identipy = args$identipy,
                  tide = args$tide, metamorpheus = args$metamorpheus,
                  msfragger = args$msfragger, maxquant = args$maxquant)

# for testing
## file_list <- list(
##   comet = "../results/test_manifest/1-First_pass/Percolator/comet_percolator_psms.tsv",
##   identipy = "../results/test_manifest/1-First_pass/Percolator/identipy_percolator_psms.tsv", maxquant = "../results/test_manifest/1-First_pass/Percolator/maxquant_percolator_psms.tsv", msfragger = "../results/test_manifest/1-First_pass/Percolator/msfragger_percolator_psms.tsv"
## )


flashlfq <- read_engine_psms(file_list[[1]], names(file_list)[1], mapping)
# Set up the first set of psms
unique_peps <- flashlfq$`Base Sequence`
#
# Add only new psms to the first set
for (n in seq_along(file_list)[-1]) {
  current <- read_engine_psms(file_list[[n]], names(file_list)[n], mapping)
  current_peps <- current$`Base Sequence`
  new_rows <- current %>% filter(!(`Base Sequence` %in% unique_peps))
  flashlfq <- bind_rows(flashlfq, new_rows)
  unique_peps <- c(
    unique_peps,
    current_peps[!(current_peps %in% unique_peps)]
  )
  rm(current, current_peps, new_rows)
}

flashlfq <- flashlfq %>%
  filter(!(is.na(`Peptide Monoisotopic Mass`))) %>%
  mutate(`Base Sequence` = unlist(lapply(`Base Sequence`, gsub, pattern = "X", replacement = "")))

write_delim(flashlfq, args$output, delim = "\t")
