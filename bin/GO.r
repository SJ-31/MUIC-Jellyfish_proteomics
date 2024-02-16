#' Helper functions for analyzing
#' and semantic similarity
#'

library(paletteer)
ONTOLOGIES <- c("MF", "BP", "CC")
library(GO.db)
library(reticulate)
library(glue)
library(GOSemSim)
library(vegan)
library(clusterProfiler)
library(tidyverse)


#' Related GO terms based on ontology structure
#'
#' @description
#' Returns GO term ids related to query term id "term",
#' which include secondary ids and all offspring (children and their children)
#' of "term"
goRelated <- function(term) {
  if (!is.null(query <- GO.db::GOTERM[[term]])) {
    ontology <- query@Ontology
    secondary <- query@Secondary
    offspring <- switch(ontology,
                        CC = { GOCCOFFSPRING[[term]] },
                        BP = { GOBPOFFSPRING[[term]] },
                        MF = { GOMFOFFSPRING[[term]] })
    return(c(secondary, offspring))
  }
  return(NULL)
}

#' Converts a matrix into a tibble, moving the row names into a column
m2Tb <- function(matrix, first_col) {
  return(as.data.frame(matrix) %>%
           tibble::rownames_to_column(var = first_col) %>%
           as_tibble())
}

#' Saves both 3d or 2d plots
mySaveFig <- function(fig, filename) {
  if (any(grepl("ggplot", class(fig)))) {
    ggsave(filename = filename, plot = fig)
  } else {
    htmlwidgets::saveWidget(fig, filename, selfcontained = TRUE)
  }
}

cosineDissimilarity <- function(tb, name_col) {
  mat <- as.matrix(tb %>% dplyr::select(., -!!name_col))
  sim <- mat / sqrt(rowSums(mat * mat))
  sim <- sim %*% t(sim)
  mat <- as.dist(1 - sim)
  return(mat)
}

#' Compute euclidean distance matrix
#'
#' @description
#' Wrapper function for R's dist function, returns a tibble
euclideanDistance <- function(tb, name_col) {
  return(
    dplyr::select(tb, where(is.numeric)) %>%
      dist(method = "euclidean") # %>%
    # as.matrix() %>%
    # as_tibble() %>%
    # mutate(!!name_col := tb[[name_col]], .before = where(is.numeric))
  )
}

#' Equalize the number of proteins between the sample and taxa being compared
#' 
#' @description
#' When a given taxon has more proteins than the sample, the proteins
#' are sorted by length and selection is performed at regular intervals
#' so that the taxon tibble has the same dimensions as the sample tb
adjustProteinNumbers <- function(sample_tb, taxa_tb_list) {
  sample_num <- nrow(sample_tb)
  nested <- taxa_tb_list %>% group_by(taxon) %>% nest()
  nested$data <- nested$data %>% lapply(\(x) {
    nrows <- nrow(x)
    if (nrows <= sample_num) {
      return(x)
    }
    step <- nrows / sample_num
    wanted_indices <- round(seq(1, nrows, by = step))
    data <- x %>%
      arrange(by = length) %>%
      dplyr::slice(wanted_indices)
    return(data)
  })
  return(unnest(nested))
}


t2Df <- function(tibble, row_col) {
  # Converts a tibble into a dataframe, with column "row_col" becoming the rownames
  converted <- tibble %>%
    dplyr::select(-!!as.symbol(row_col)) %>%
    as.data.frame() %>%
    `row.names<-`(tibble[[row_col]])
  return(converted)
}


#' Label row indices that contain a string query
#' in the specified column
#'
#' @description
#' If the query isn't found in the target column,
#' then NA is returned instead for that entry
markMatch <- function(tb, target_col, query, label) {
  purrr::pluck(tb, target_col) %>%
    purrr::map_chr(\(x) if (grepl(query, x, ignore.case = TRUE)) label else NA)
}


goVector <- function(df, column, filter, go_column) {
  if (!missing(go_column)) {
    df <- df %>% dplyr::rename(GO_temp = go_column)
  }
  if (!missing(column)) {
    df <- df %>% dplyr::filter(!!as.symbol(column) == filter)
  }
  gos <- df$GO_temp %>%
    lapply(., str_split_1, pattern = ";") %>%
    unlist() %>%
    discard(is.na) %>%
    discard(!grepl("GO:", .))
  return(gos)
}

prepOrgDb <- function(path) {
  name <- basename(path)
  install.packages(path, repos = NULL)
  library(name, character.only = TRUE)
}


#'  Reduce a  list of GO terms by removing redundant terms
#'  Returns three types of results: the matrix of representative go terms,
#'  the similarity matrix and the associated scatter plot
#'  There are three for each ontology in the GO
reduceGOList <- function(go_list) {
  sims <- lapply(ONTOLOGIES, \(x) {
    rrvgo::calculateSimMatrix(go_list,
                              orgdb = db_name,
                              keytype = "GID",
                              ont = x,
                              method = "Rel"
    )
  })
  names(sims) <- ONTOLOGIES
  reduce_matrices <- lapply(sims, \(x) {
    close <- myReduceSimMatrix(
      x,
      threshold = 0.7,
      orgdb = db_name,
      keytype = "GID"
    )
    return(as_tibble(close))
  })
  names(reduce_matrices) <- ONTOLOGIES
  plots <- lapply(ONTOLOGIES, \(x) {
    scatterPlot(sims[[x]], reduce_matrices[[x]])
  }) %>% `names<-`(ONTOLOGIES)
  return(list(
    reduced_matrix = reduce_matrices,
    sim_matrix = sims,
    scatter_plots = plots
  ))
}


makeDistMatrix <- function(v1, v2, dist_func) {
  # Construct a distance matrix from vectors v1 and v2 by applying "dist_func"
  # to each pair of elements
  return(sapply(v1, \(x) {
    sapply(v2, \(y) dist_func(y, x))
  }))
}

### PCOA
pcoaWithTb <- function(distances, join_on) {
  # Perform principal coordinates analysis
  # Merge with the annotation tibble
  pcoa <- vegan::wcmdscale(distances, eig = TRUE)
  return(pcoa)
}


ontoResults <- function(ontologizer_dir) {
  # Convenience function for aggregating statistically significant ontologizer results
  # There will be two different tibbles, containing GO terms known only
  #    from...
  #    - de novo and transcriptome peptides
  #    - with detected PTMs
  onto_files <- list.files(ontologizer_dir) %>%
    keep(., grepl("ontologizer-", .)) %>%
    lapply(., \(x) {
      glue("{ontologizer_dir}{x}")
    }) %>%
    unlist()
  names(onto_files) <- onto_files %>%
    lapply(., str_match,
           pattern = ".*-(.*)\\.txt"
    ) %>%
    lapply(., \(x) x[, 2]) %>%
    unlist()
  onto_files <- lapply(onto_files, readr::read_tsv)
  onto_files$gos_from_downloads <- onto_files$from_downloaded_db %>%
    filter(!is.trivial) %>%
    goVector(go_column = "ID")
  onto_files$gos_from_open <- onto_files$id_with_open %>%
    filter(!is.trivial) %>%
    goVector(go_column = "ID")
  return(onto_files)
}

# normalizeGoCounts <- function(nested)


goDataGlobal <- function(uniprot_data_dir, sample_data, sample_name, onto_path) {
  # Load sample data & ontologizer results
  sample_tb <- readr::read_tsv(sample_data) %>%
    filter(!is.na(GO_IDs)) %>%
    mutate(ProteinId = paste0(ProteinId, "-SAMPLE"))
  ontologizer <- ontoResults(onto_path)

  # Load Uniprot data into a list of tibbles
  file_list <- list.files(uniprot_data_dir, "*_reviewed.tsv",
                          full.names = TRUE
  )
  all_tbs <- file_list %>%
    lapply(\(x) {
      read_tsv(x) %>%
        mutate(
          `Gene Ontology IDs` = unlist(lapply(
            `Gene Ontology IDs`, gsub,
            pattern = " ", replacement = "", .
          )),
          taxon = gsub(".*/(.*)_reviewed\\.tsv", "\\1", x)
        ) %>%
        dplyr::rename(length = Length)
    }) %>%
    `names<-`(lapply(file_list, gsub,
                     pattern = ".*/(.*)\\.tsv", replacement = "\\1"
    ))


  # Extract all GO terms from the list of tibbles
  # Load into new combined tibble and merge with sample GO terms
  # Load sample GO terms into separate tibble
  go_tb_all <- names(all_tbs) %>%
    lapply(., \(x) {
      tb <- all_tbs[[x]] %>% filter(!is.na(GO_IDs))
      gos <- goVector(df = tb, go_column = "GO_IDs")
      return(tibble(GO_IDs = gos, taxon = gsub("_reviewed", "", x)))
    }) %>%
    bind_rows()

  sample_gos <- goVector(sample_tb, go_column = "GO_IDs") %>% unique()

  go_tb_all <- bind_rows(go_tb_all, tibble(
    GO_IDs = sample_gos,
    taxon = sample_name
  ))
  # Taxa that have a higher (normalized) frequency of a given GO ID will be
  # assigned that id
  # TODO: Is this really the best way to do it? Who could you ask?
  nested <- go_tb_all %>%
    group_by(GO_IDs) %>%
    nest()
  taxa_counts <- go_tb_all$taxon %>% table()
  go_tb_all <- lapply(seq(nrow(nested)), \(x) {
    go_id <- nested[x,]$GO_IDs
    cur_nest <- nested[x,]$data[[1]] %>% table()
    normalized <- cur_nest / taxa_counts[names(taxa_counts) %in% names(cur_nest)]
    most_frequent <- normalized %>%
      as_tibble() %>%
      arrange(desc(n)) %>%
      dplyr::slice(1)
    return(tibble(GO_IDs = go_id, taxon = most_frequent$taxon))
  }) %>% bind_rows()

  # Load four vectors of GO terms:
  # 1. Sample, 2. All go terms, 3. GO terms from
  # eggnog proteins 4. GO terms in interpro proteins

  all_gos <- c(go_tb_all$GO_IDs, sample_gos) %>% unique()
  interpro_gos <- goVector(sample_tb, go_column = "GO_IDs", "Anno_method", "interpro")
  eggnog_gos <- goVector(sample_tb, go_column = "GO_IDs", "Anno_method", "eggNOG")

  go_tb_sample <- tibble(
    GO_IDs = sample_gos,
    sig_downloaded_db = sample_gos %in% ontologizer$gos_from_downloads,
    sig_id_w_open = sample_gos %in% ontologizer$id_with_open
  )
  prot_tb_sample <- dplyr::select(sample_tb, c("ProteinId", "GO_IDs"))
  prot_tb_all <- all_tbs %>%
    bind_rows() %>%
    dplyr::rename(ProteinId = Entry) %>%
    bind_rows(dplyr::mutate(sample_tb, taxon = sample_name)) %>%
    dplyr::filter(!is.na(GO_IDs))
  data <- list()
  data$sample_tb <- sample_tb
  data$ontologizer <- ontologizer
  data$all_tbs <- all_tbs
  data$go_tb$all <- go_tb_all
  data$go_tb$sample <- go_tb_sample
  data$go_vec$sample <- sample_gos
  data$go_vec$all <- all_gos
  data$go_vec$interpro <- interpro_gos
  data$go_vec$eggnog <- eggnog_gos
  data$protein$all <- adjustProteinNumbers(prot_tb_sample, prot_tb_all)
  data$protein$sample <- prot_tb_sample
  data$prot_go_map$all <- prot_tb_all$GO_IDs %>%
    lapply(., str_split_1, pattern = ";") %>%
    `names<-`(prot_tb_all$ProteinId)
  data$prot_go_map$sample <- data$prot_go_map$all %>%
    purrr::keep_at(\(x) x %in% prot_tb_sample$ProteinId)
  return(data)
}


goQueryProt <- function(query_gos, tb, prot_go_map) {
  # Queries a tibble to return only entries whose GO terms contain every GO term
  # in "query_gos" i.e. returns entries iff "query_gos" is a proper subset
  # of an entry's go terms
  wanted_ids <- prot_go_map %>%
    keep(\(x) all(query_gos %in% x)) %>%
    names()
  return(dplyr::filter(tb, ProteinId %in% wanted_ids))
}


# TODO: Ranked gene set enrichment analysis
prepSorted <- function(df, col_spec, type) {
  select_vals <- df %>%
    dplyr::select(contains(glue("{col_spec}|ProteinId"))) %>%
    rowwise(ProteinId)
  if (type == "mean") {
    calc <- select_vals %>% reframe(mean = mean(c_across(where(is.numeric))))
  } else {
    calc <- reframe(median = median(c_across(where(is.numeric))))
  }

  vec <- unlist(calc[, 2], FALSE) %>% `names<-`(unlist(calc[, 1], FALSE))
  return(sort(vec, decreasing = TRUE))
}


# Analyses list
#   Representative GO terms
