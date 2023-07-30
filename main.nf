include { search } from "./subworkflows/shotgun.nf"
include { assemble } from "./subworkflows/rna_seq.nf"

workflow rnaseq {
    assemble()
}

workflow identify {
    search()
}

