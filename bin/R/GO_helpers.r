#' Helper functions for analyzing GO terms
#' and semantic similarity
#'

library("paletteer")
library("tidyverse")
ONTOLOGIES <- c("MF", "BP", "CC")
library("GO.db")
library("reticulate")
library("fgsea")
library("GOSemSim")
library("vegan")
library("clusterProfiler")

#' Unique values of a named vector (based on its name)
#'
uniqueNames <- function(vec) {
  return(vec[!duplicated(names(vec))])
}


#' Related GO terms based on ontology structure
#'
#' @description
#' Returns GO term ids that related to query term id "term",
#' which includes secondary ids and all offspring (children and their children)
#' of "term"
goOffspring <- function(term) {
  if (!is.null(query <- GO.db::GOTERM[[term]])) {
    ontology <- query@Ontology
    secondary <- query@Secondary
    offspring <- switch(ontology,
      CC = {
        GOCCOFFSPRING[[term]]
      },
      BP = {
        GOBPOFFSPRING[[term]]
      },
      MF = {
        GOMFOFFSPRING[[term]]
      }
    )
    return(c(secondary, offspring))
  }
  return(NULL)
}

countOffspring <- function(term) {
  return(length(goOffspring(term)))
}


#' Get semantic value of  a GO term (from Wang's method)
#'
#' @description
#' A wrapper function around the non-exported getSV function from GOSemSim
getSV <- function(term) {
  if (!exists(".GOSemSimEnv")) {
    init <- getFromNamespace(".initial", "GOSemSim")
    init()
  }
  t <- GOTERM[[term]]
  if (is.null(t)) {
    return(NA)
  }
  ont <- t@Ontology
  rel_df <- base::get("gotbl", envir = .GOSemSimEnv)
  sv <- getFromNamespace("getSV", "GOSemSim")(term, ont, rel_df)
  return(sum(sv))
}


#' Returns TRUE if term a is offspring of term b
#'
isGoOffspring <- function(a, b) {
  offspring <- goOffspring(b)
  if (is.null(offspring)) {
    return(FALSE)
  }
  return(a %in% offspring)
}

#' Convert dataframe to tibble, moving row names into a named column
df2Tb <- function(df, first_col) {
  return(rownames_to_column(df, var = first_col) %>% as_tibble())
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
    extras <- gsub(".html", "_files", filename)
    unlink(extras, recursive = TRUE)
  }
}


#' Equalize the number of proteins between the sample and taxa being compared
#'
#' @description
#' When a given taxon has more proteins than the sample, the proteins
#' are sorted by length and selection is performed at regular intervals
#' so that the taxon tibble has the same dimensions as the sample tb
adjustProteinNumbers <- function(target_num, tb) {
  nested <- tb %>%
    group_by(Taxon) %>%
    nest()
  picked <- nested %>%
    apply(
      1,
      \(x) {
        data <- x$data
        nrows <- nrow(data)
        if (nrows <= target_num) {
          return(data)
        }
        data <- data[sample(target_num), ]
        return(slice(data, 1:target_num))
      }
    ) %>%
    bind_rows()
  return(picked)
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


#' Retrieve a vector of GO terms from a tibble where terms for a given
#' protein have been combined by ";"
#'
#' @description
#' @param go_column the column of `tb` containing the terms
#' @param filter optional filter criteria to use with `column`
goVector <- function(tb, column, filter, go_column, unique = FALSE) {
  if (any(is.na(tb[[go_column]]))) {
    tb <- tb %>% dplyr::filter(!is.na(!!as.symbol(go_column)))
  }
  if (!missing(go_column)) {
    tb <- tb %>% dplyr::rename(GO_temp = go_column)
  }
  if (!missing(column)) {
    tb <- tb %>% dplyr::filter(!!as.symbol(column) == filter)
  }
  gos <- tb$GO_temp %>%
    lapply(., str_split_1, pattern = ";") %>%
    unlist() %>%
    discard(is.na) %>%
    discard(!grepl("GO:", .))
  if (unique) {
    return(unique(gos))
  }
  return(gos)
}

prepOrgDb <- function(path) {
  name <- base::basename(path)
  if (!name %in% installed.packages()[, 1]) {
    try(install.packages(path, repos = NULL))
  }
  library(name, character.only = TRUE)
  name
}

#'  Reduce a  list of GO terms by removing redundant terms
#' Uses the specified semantic similarity method to generate distance
#' matrix, then performs hierarchical clustering
#'  Returns three types of results: the matrix of representative go terms,
#'  the similarity matrix and the associated scatter plot
#'  There are three for each ontology in the GO
#' Without scores, it relies on term uniqueness
#' @param go_list A named vector (names are GO terms) of scores
reduceGOList <- function(go_list) {
  sims <- lapply(ONTOLOGIES, \(x) {
    rrvgo::calculateSimMatrix(names(go_list),
      orgdb = DB_NAME,
      keytype = "GID",
      ont = x,
      method = "Wang" # Because the IC-based methods are based on the
      # specific corpus of GO terms, this can introduce bias
      # Wang's graph-based method does not have this limitation
    )
  })
  names(sims) <- ONTOLOGIES
  # The scores generated here are from hiearchichal clustering
  reduce_matrices <- lapply(sims, \(x) {
    close <- myReduceSimMatrix(
      x,
      threshold = 0.9,
      orgdb = DB_NAME,
      scores = go_list,
      keytype = "GID"
    )
    return(as_tibble(close))
  })
  names(reduce_matrices) <- ONTOLOGIES
  plots <- lapply(ONTOLOGIES, \(x) {
    scatterPlot(sims[[x]], reduce_matrices[[x]])
    # Applies PCoA to the similarity matrix
  }) %>% `names<-`(ONTOLOGIES)
  return(list(
    reduced_matrix = reduce_matrices,
    sim_matrix = sims,
    scatter_plots = plots
  ))
}

termOntology <- function(term) {
  ontology <- ifelse(is.null(GOTERM[[term]]),
    "NONE", GOTERM[[term]]@Ontology
  )
  return(ontology)
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
      glue("{ontologizer_dir}/{x}")
    }) %>%
    unlist()
  names(onto_files) <- onto_files %>%
    lapply(., str_match,
      pattern = ".*-(.*)\\.txt"
    ) %>%
    lapply(., \(x) x[, 2]) %>%
    unlist()
  onto_files <- lapply(onto_files, readr::read_tsv)
  getNamedGO <- function(tb) {
    filtered <- dplyr::filter(tb, !is.trivial)
    filtered %>%
      purrr::pluck("p.adjusted") %>%
      `names<-`(filtered$ID)
  }
  # Named vector of GO terms where values are adjusted p_values
  onto_files$unknown_to_db_GO <- getNamedGO(onto_files$unknown_to_db)
  onto_files$id_with_open_GO <- getNamedGO(onto_files$id_with_open)
  return(onto_files)
}


#' Create a tibble containing information about specific GO terms
#'
goInfoTb <- function(go_vector) {
  assertArg(go_vector, \(x) is.atomic(x))
  tb <- lapply(go_vector, \(x) {
    row <- tibble(
      GO_IDs = x, term = NA,
      definition = NA, ontology = NA,
    )
    find_info <- GOTERM[[x]]
    if (is.null(find_info)) {
      return(row)
    }
    row$term <- find_info@Term
    row$ontology <- find_info@Ontology
    row$definition <- find_info@Definition
    return(row)
  }) %>% bind_rows()
  return(tb)
}

getUniprotData <- function(uniprot_tsv_path) {
  tb <- read_tsv(uniprot_tsv_path) %>%
    filter(!is.na(GO_IDs)) %>%
    rename(ProteinId = Entry)
  go_tb_all <- tb %>%
    group_by(Taxon) %>%
    nest() %>%
    apply(
      1,
      \(x) {
        taxon <- x$Taxon
        cur <- x$data
        return(tibble(
          GO_IDs = goVector(cur, go_column = "GO_IDs"),
          Taxon = taxon
        ))
      }
    ) %>%
    bind_rows()
  all_gos <- go_tb_all$GO_IDs %>% unique()
  prot_go_map <- tb$GO_IDs %>%
    lapply(., str_split_1, pattern = ";") %>%
    `names<-`(tb$ProteinId)
  return(list(
    map = prot_go_map,
    go_vec = all_gos,
    prot_tb = ungroup(tb),
    go_tb = go_tb_all
  ))
}

goData <- function(sample_path, onto_path, uniprot_tsv, sample_name) {
  # Load sample data & ontologizer results
  sample_tb <- readr::read_tsv(sample_path) %>%
    filter(!is.na(GO_IDs)) %>%
    distinct(ProteinId, .keep_all = TRUE)
  sample_gos <- goVector(sample_tb, go_column = "GO_IDs") %>% unique()

  if (!missing(uniprot_tsv) && !is.null(uniprot_tsv)) {
    # Load Uniprot data into a list of tibbles
    unip <- getUniprotData(uniprot_tsv)
    # Extract all GO terms from the list of tibbles
    # Load into new combined tibble and merge with sample GO terms
    # Load sample GO terms into separate tibble
    go_tb_all <- bind_rows(unip$go_tb, tibble(
      GO_IDs = sample_gos,
      Taxon = sample_name
    ))
    # Taxa that have a higher (normalized) frequency of a given GO ID will be
    # assigned that id
    nested <- go_tb_all %>%
      group_by(GO_IDs) %>%
      nest()
    taxa_counts <- go_tb_all$Taxon %>% table()
    go_tb_all <- lapply(seq(nrow(nested)), \(x) {
      go_id <- nested[x, ]$GO_IDs
      cur_nest <- nested[x, ]$data[[1]] %>% table()
      normalized <- cur_nest / taxa_counts[names(taxa_counts) %in% names(cur_nest)]
      if (dim(normalized) == 0) {
        return(tibble())
      }
      most_frequent <- normalized %>%
        as_tibble() %>%
        arrange(desc(n)) %>%
        dplyr::slice(1)
      return(tibble(GO_IDs = go_id, Taxon = most_frequent$Taxon))
    }) %>% bind_rows()

    all_gos <- c(go_tb_all$GO_IDs, sample_gos) %>% unique()
  }

  # Load four vectors of GO terms:
  # 1. Sample, 2. All go terms, 3. GO terms from
  # eggnog proteins 4. GO terms in interpro proteins

  interpro_gos <- goVector(sample_tb, go_column = "GO_IDs", "inferred_by", "interpro")
  eggnog_gos <- goVector(sample_tb, go_column = "GO_IDs", "inferred_by", "eggNOG")

  prot_tb_sample <- dplyr::select(sample_tb, c("ProteinId", "GO_IDs"))
  data <- list()
  data$sample_tb <- sample_tb
  data$go_vec$sample <- sample_gos
  data$go_vec$interpro <- interpro_gos
  data$go_vec$eggnog <- eggnog_gos
  data$protein$sample <- prot_tb_sample
  data$prot_go_map$sample <- prot_tb_sample$GO_IDs %>%
    lapply(., str_split_1, pattern = ";") %>%
    `names<-`(prot_tb_sample$ProteinId)

  if (!missing(onto_path)) {
    ontologizer <- ontoResults(onto_path)
    go_tb <- tibble(
      GO_IDs = sample_gos,
      sig_downloaded_db = sample_gos %in% ontologizer$gos_from_downloads,
      sig_id_w_open = sample_gos %in% ontologizer$id_with_open
    )
    data$ontologizer <- ontologizer
    data$go_tb$sample <- go_tb
  }
  if (!missing(uniprot_tsv) && !is.null(uniprot_tsv)) {
    data$protein$compare <- unip$prot_tb
    data$prot_go_map$compare <- unip$map
    data$go_vec$compare <- all_gos
    data$go_tb$compare <- go_tb_all
    data$unip_tb <- unip$all_tb
  }
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


TARGET_TERMS <- list(toxins = c(
  "GO:0090719", # Toxin activity
  "GO:0046930", # pore complex/pore-forming toxin activity
  "GO:0015288", # porin activity
  "GO:0031640", # Killing cells of another organism
  "GO:0015473" # Fimbrial usher porin activity
))
TARGET_TERMS <- purrr::map(TARGET_TERMS, \(x) {
  purrr::map(x, goOffspring) %>%
    unlist() %>%
    unique() %>%
    discard(is.na)
})

getToxinProteins <- function(prot2go_map) {
  prot2go_map %>%
    keep(\(x) any(x %in% TARGET_TERMS$toxins))
}

# A list of important, higher-level GO categories (somewhat arbitrary)
# Can have up to 20 to be visualized (due to restrictions with
# color palettes)
GO_CATEGORIES <- list(
  venom_component = c(
    "GO:0090719", "GO:0046930", # 1
    "GO:0031640", "GO:0015473"
  ),
  small_molecule_binding = "GO:0036094", # 2
  translation = "GO:0006412", # 3
  transport = "GO:0006810", # 4
  cell_projection = "GO:0042995", # 5
  cytoskeleton = "GO:0005856", # 6
  membrane = "GO:0016020", # 7
  catalytic_activity = "GO:0003824", # 8
  organelle = "GO:0043226" # 9
)
GO_CATEGORIES <- lapply(GO_CATEGORIES, \(x) {
  offspring <- lapply(x, goOffspring) %>% unlist()
  return(c(x, offspring))
})

headerFreqs <- function(header_vec) {
  header_vec %>%
    lapply(., str_split_1, pattern = " ") %>%
    unlist() %>%
    table() %>%
    sort(decreasing = TRUE) %>%
    as_tibble() %>%
    rename(c("word" = ".", "count" = "n")) %>%
    filter(
      !grepl("=", word),
      nchar(word) > 2,
      grepl("[A-Za-z]+", word),
      !grepl("\\.|,", word)
    )
  # Remove database classifications
  # Get rid of tiny words
}




GO_DAG <- NULL
GO_SLIM_DAG <- NULL
getGoSlim <- function(go_terms, go_path, go_slim_path, which = NULL) {
  op <- reticulate::import("goatools.obo_parser")
  ms <- reticulate::import("goatools.mapslim")
  if (is.null(GO_DAG)) {
    GO_DAG <<- op$GODag(go_path)
    GO_SLIM_DAG <<- op$GODag(go_slim_path)
  }
  getOneSlim <- function(term) {
    py$slims <- tryCatch(
      expr = ms$mapslim(term, go_dag = GO_DAG, goslim_dag = GO_SLIM_DAG),
      error = function(cnd) NULL
    )
    if (is.null(py$slims) || length(py$slims) == 0) {
      return(list(direct = NA, all = NA))
    }
    reticulate::py_run_string("direct = list(slims[0])")
    reticulate::py_run_string("all = list(slims[1])")
    result <- list(direct = py$direct, all = py$all)
    if (!is.null(which) && which == "direct") {
      return(result$direct)
    } else if (!is.null(which) && which == "all") {
      return(result$all)
    } else {
      return(result)
    }
  }
  slims <- lapply(go_terms, getOneSlim) %>% `names<-`(go_terms)
  return(slims)
}

#' Split a string of GO terms joined by ";" and return their slims joined together
slimsFromGoString <- function(term_str, go_path, go_slim_path) {
  if (is.na(term_str)) {
    return(NA)
  }
  slims <- str_split_1(term_str, ";") %>%
    lapply(
      ., \(g) getGoSlim(g, go_path, go_slim_path, "all")
    ) %>%
    unlist() %>%
    unique() %>%
    discard(is.na) %>%
    paste0(collapse = ";")
  if (slims == "") {
    return(NA)
  }
  return(slims)
}

#' Find the terms from a vector of GO Ids, grouping them by sub-ontology
#' @return a list of three GO term vectors, one for each sub-ontology
idsIntoOntology <- function(id_vector, target = "Term", collapse = TRUE) {
  assertArg(id_vector, \(x) is.atomic(x) || is.null(x))
  empty <- list("CC" = "", "BP" = "", "MF" = "")
  if (is.null(id_vector)) {
    return(empty)
  }
  getTerm <- function(id) {
    list <- empty
    term <- GOTERM[[id]]
    if (!is.null(term)) {
      attrs <- attributes(term)
      list[[term@Ontology]] <- attrs[[target]]
    }
    return(list)
  }
  reduced <- lapply(id_vector, getTerm) %>%
    purrr::reduce(., mergeLists, .init = empty)
  if (collapse) {
    reduced <- purrr::map(reduced, \(x) {
      x %>%
        discard(., x == "") %>%
        paste0(., collapse = ";")
    })
  }
  return(reduced)
}

# Map vector of Protein ids their groups,
mapUnique <- function(ids, map) {
  map_chr(ids, \(x) map[[x]]) %>%
    unique() %>%
    discard(is.na) %>%
    discard(\(x) x == "U")
}

getOrganism <- function(tb) {
  organismFromHeader <- function(row) {
    # Parse the organism from an NCBI or uniprot ID
    organism <- row["organism"]
    if (is.na(organism)) {
      header <- row["header"]
      if (grepl("OS=", header)) {
        return(str_extract(header, "OS=([a-zA-Z]* [a-zA-Z]*) ", group = 1))
      }
      return(str_extract(header, ".*\\[([A-Z].*)\\]", group = 1))
    }
    return(organism)
  }
  tb |> mutate(organism = apply(tb, 1, organismFromHeader))
}
