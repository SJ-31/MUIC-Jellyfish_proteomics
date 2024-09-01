library("tidyverse")
library("glue")

b_header <- "queryId,subjectId,sAlignStart,sAlignEnd,alignLen,bitscore,evalue,pident,nident,nmismatch,ngaps"

if (str_detect(getwd(), "Bio_SDD")) {
  wd <- "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
} else {
  wd <- "/home/shannc/workflow"
}
if (!exists("SOURCED")) {
  source(glue("{wd}/bin/R/analysis/all_analyses.r"))
  SOURCED <- TRUE
}
data <- read_tsv(M$data_w_cat_path)

qdir <- glue("{wd}/docs/figures/query_seqs")
bdir <- glue("{wd}/docs/figures/query_seqs/blast")

id_from_header <- function(header) {
  if (str_detect(header, "\\|")) {
    return(str_extract(header, ".*\\|([A-Z+0-9]+)\\|", group = 1))
  }
  header
}

join_with_named <- function(tb) {
  tb <- tb |> mutate(
    subjectId = map_chr(subjectId, id_from_header),
    row_id = seq_len(nrow(tb))
  )
  joined <- list()
  for (col in c("NCBI_ID", "UniProtKB_ID", "header")) {
    D <- data |> mutate(join = !!as.symbol(col))
    cur <- inner_join(tb, D, by = join_by(x$subjectId == y$join))
    cur$joined <- TRUE
    joined[[col]] <- cur
    tb <- tb |> filter(!row_id %in% cur$row_id)
  }
  joined <- joined |>
    discard(\(x) nrow(x) == 0) |>
    bind_rows()
  bind_rows(joined, mutate(tb, joined = FALSE))
}

write_named <- function(file) {
  name <- str_remove(file, ".*/") |> str_remove("_blast\\.csv")
  n <- read_csv(file, col_names = header) |>
    join_with_named() |>
    arrange(desc(pident)) |>
    select(-c(ProteinGroupId, GroupSB, sb_rep))
  n |> write_tsv(glue("{qdir}/named/{name}_named.tsv"))
  n
}

results_files <- list.files(bdir, pattern = ".*_blast.csv", full.names = TRUE)
results_files
header <- str_split_1(b_header, ",")
nema <- write_named(results_files[7])
nowa <- write_named(results_files[8])
minicollagen <- write_named(results_files[6])
