library(ggkegg)
library(KEGGREST)

#' Node vectors
#'
#' @description
#' Return all nodes from a ggkegg igraph object
getNodes <- function(g) {
  if (is.null(g)) {
    return(NA)
  } else {
    return(
      g %>%
        activate(nodes) %>%
        as_tibble() %>%
        pluck("name") %>%
        flattenJoined(., " ")
    )
  }
}
