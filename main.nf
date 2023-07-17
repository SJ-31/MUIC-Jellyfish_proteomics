include { database_search } from "./subworkflows/shotgun.nf"
include { de_novo } from "./subworkflows/shotgun.nf"
include { assemble } from "./subworkflows/rna_seq.nf"

workflow rnaseq {
    assemble()
}

workflow search {
    database_search()
}

workflow de_novo {
    de_novo()
}
