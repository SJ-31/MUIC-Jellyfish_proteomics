library(prozor)
args <- commandArgs(trailingOnly = TRUE)
decoy <- createDecoyDB(dbs = args[1], revLab = "rev_")
writeFasta(decoy, file = "./all_decoy.fasta")
