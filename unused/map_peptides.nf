process MAP_PEPTIDES {
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_results)
    val(outdir)
    //

    output:
    path("peptides2proteinId.tsv")
    //

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(glue)

    source(glue("${params.bin}/R/helpers.r"))
    combined <- read_tsv("${combined_results}")
    mapping <- combined %>%
        filter(!is.na(MatchedPeptideIds)) %>%
        select(c(MatchedPeptideIds, ProteinId)) %>%
        tibbleDuplicateAt( "MatchedPeptideIds", ";")

    write_tsv(mapping, "peptides2proteinId.tsv")
    """
    //
}
