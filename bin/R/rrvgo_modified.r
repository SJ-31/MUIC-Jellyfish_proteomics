library(rrvgo)

##   Sayols, S (2023). rrvgo: a Bioconductor package for interpreting
##   lists of Gene Ontology terms. microPublication Biology.
##   10.17912/micropub.biology.000811
##
##   "getGoSize" function has been modified to work with custom org.db packages
##   (those without ENTREZID columns)


myGetGoSize <- function(terms, orgdb, keytype, children) {
  if (all(is(orgdb) != "OrgDb")) {
    orgdb <- rrvgo:::loadOrgdb(orgdb)
  }
  go <- AnnotationDbi::select(orgdb, keytype = if (children) {
    "GOALL"
  } else {
    "GO"
  }, keys = terms, columns = keytype)
  counts <- tapply(go$GOALL, go$GO, function(x) length(unique(x)))
  empty <- terms[!(terms %in% names(counts))]
  nocounts <- setNames(rep(0, length(empty)), empty)
  c(counts, nocounts)
}

myReduceSimMatrix <- function(
    simMatrix, scores = c("uniqueness", "size"), threshold = 0.7,
    orgdb, keytype = "ENTREZID", children = TRUE) {
  if (is(scores, "character")) {
    scores <- match.arg(scores)
  } else {
    stopifnot(`Scores vector does not contain all terms in the similarity matrix` = all(rownames(simMatrix) %in%
      names(scores)))
  }
  cluster <- cutree(hclust(as.dist(1 - simMatrix)), h = threshold)
  sizes <- myGetGoSize(rownames(simMatrix), orgdb, keytype, children)
  termUniq <- rrvgo:::getTermUniq(simMatrix, cluster)
  if (is(scores, "character")) {
    scores <- switch(scores,
      uniqueness = {
        message("No scores provided. Falling back to term's uniqueness")
        termUniq
      },
      size = {
        message("No scores provided. Falling back to term's GO size")
        sizes
      }
    )
  }
  scores <- scores[match(rownames(simMatrix), names(scores))]
  o <- rev(order(scores, termUniq, na.last = FALSE))
  scores <- scores[o]
  simMatrix <- simMatrix[o, o]
  cluster <- cluster[match(rownames(simMatrix), names(cluster))]
  clusterRep <- tapply(
    rownames(simMatrix), cluster,
    function(x) x[which.max(scores[x])]
  )
  data.frame(
    go = rownames(simMatrix),
    cluster = cluster,
    parent = clusterRep[cluster],
    score = scores[match(rownames(simMatrix), names(scores))],
    size = sizes[match(rownames(simMatrix), names(sizes))],
    term = rrvgo:::getGoTerm(rownames(simMatrix)),
    parentTerm = rrvgo:::getGoTerm(clusterRep[cluster]),
    termUniqueness = rrvgo:::getTermUniq(simMatrix),
    termUniquenessWithinCluster = rrvgo:::getTermUniq(
      simMatrix,
      cluster
    ),
    termDispensability = rrvgo:::getTermDisp(simMatrix, cluster, clusterRep)
  )
}
