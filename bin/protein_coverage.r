library(tidyverse)
library(msa)

prot_table <- "../results/ref/all_normal_mapping.tsv"
seq_mapping <- "../results/ref/all_normal_mapping.tsv"

prot_df <- read.delim(prot_table, sep = "\t")
mapping <- read.delim(seq_mapping, sep = "\t")

myprot <-  "AFKDJFWFSJVSPPFSAVSESFFSEFWWSSSFEQVXCFEQFDSFSFEAFSFEREFSAFDSWWJJSFS"
mypeps <- c("WFSJVSP", "FWWSSSFEQVXC", "FSFEREFSAFDSWWJJSF")


protein <- gsub("[^A-Z]+", "", myprot)
mask <- rep(0, nchar(protein))

all_locs <- lapply(mypeps, function(x) {
  return(str_locate_all(protein, x)[[1]])
})

sapply(1:dim(locs)[1], function(x) {
  start <- locs[x, "start"] %>% unlist(use.names = FALSE)
  mask[unname(locs[x, "start"]):unname(locs[x, "end"])] <<- 1
})

mask
