library("tidyverse")
library("rhdf5")

orgdb <- prepOrgDb(args$orgdb_path)

ALL_SEM_DATA <- list(
  BP = GOSemSim::godata(orgdb, ont = "BP"),
  MF = GOSemSim::godata(orgdb, ont = "MF"),
  CC = GOSemSim::godata(orgdb, ont = "CC")
)

writeSim <- function(
    hdf5file, sim = ? Matrix(),
    ontology = ? Character()) {
  if (!file.exists(hdf5file)) {
    h5createFile(hdf5file)
    h5createGroup(hdf5file, "matrices")
    h5createGroup(hdf5file, "names")
  }
  h5write(sim, hdf5file, glue("matrices/{ontology}"))
  h5write(rownames(sim), hdf5file, glue("names/{ontology}"))
}

has_go <- data |> filter(!is.na(GO_IDs))
for (o in names(ALL_SEM_DATA)) {
  sim <- GOSemSim::mgeneSim(
    as.list(has_go$ProteinId),
    ALL_SEM_DATA[[o]]
  )
  writeSim(SEM_MATRIX_FILE, sim, o)
}
