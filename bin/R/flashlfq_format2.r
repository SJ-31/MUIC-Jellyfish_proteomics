library(tidyverse)
library(optparse)


main <- function(args) {
  flashlfq_header <- c(
    "File Name", "Scan Retention Time", "Precursor Charge",
    "Base Sequence", "Full Sequence", "Peptide Monoisotopic Mass",
    "Protein Accession"
  )
  old_names <- c(
    "file", "retentionTime", "precursorCharge",
    "base_peptide", "peptide", "mw", "protein"
  )

  files <- paste0(args$path, "/", list.files(args$path))

  all_engines <- lapply(files, function(x) {
    current <- read.delim(x, sep = "\t") %>%
      mutate(retentionTime = retentionTime / 60) %>%
      select(all_of(old_names)) %>%
      rename_with(~flashlfq_header, all_of(old_names)) %>%
      as_tibble()
  }) %>% bind_rows()

  all_engines <- all_engines %>%
    filter(!(is.na(`Peptide Monoisotopic Mass`))) %>%
    mutate(
      `Base Sequence` =
        unlist(lapply(`Base Sequence`, gsub, pattern = "X", replacement = "")),
      `Full Sequence` =
        unlist(lapply(`Full Sequence`, gsub,
                      pattern = "-\\.|\\.-",
                      replacement = ""
        ))
    ) %>%
    filter(!is.na(`Protein Accession`))
  grouped <- all_engines %>%
    group_by(`Base Sequence`, `File Name`, `Full Sequence`) %>%
    mutate(across(c("Scan Retention Time", "Protein Accession"), \(x) {
      paste0(x, collapse = ";")
    })) %>%
    distinct() %>%
    ungroup() %>%
    mutate(
      `Scan Retention Time` =
        unlist(lapply(`Scan Retention Time`, \(x) {
          unique(str_split_1(x, ";")) %>%
            as.double() %>%
            median()
        })),
      `Protein Accession` =
        unlist(lapply(`Protein Accession`, \(x) {
          unique(str_split_1(x, ";")) %>% paste0(collapse = ";")
        }))
    )
  return(grouped)
}

if (sys.nframe() == 0) { # Won't run if the script is being sourced
  parser <- OptionParser()
  parser <- add_option(parser, c("-p", "--path"),
                       type = "character",
                       help = "Path to formatted files"
  )
  parser <- add_option(parser, c("-o", "--output"),
                       type = "character",
                       help = "Output file name"
  )
  args <- parse_args(parser)
  o <- main(args)
  write_delim(o, args$output, delim = "\t")
}
