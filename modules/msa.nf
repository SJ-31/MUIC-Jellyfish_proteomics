process MSA {
    publishDir "$outdir", mode: "copy"
    publishDir "params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("aligned.fasta"), emit: alignment
    path("to_align.fasta")
    path("*.log")
    //

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(seqinr)
    tb <- read_tsv("${combined_tsv}")
    tb <- tb %>% filter(pcoverage_align > 0.4)
    write.fasta(as.list(tb\$seq), tb\$ProteinId, "to_align.fasta", as.string = TRUE)
    system2("mafft", c("to_align.fasta", ">", "aligned.fasta"))
    system2("cp", c(".command.err", "mafft_align.log"))
    """
    //
}
