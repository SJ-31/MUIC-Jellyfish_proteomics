include { search } from "./subworkflows/shotgun.nf"
include { assemble } from "./subworkflows/rna_seq.nf"
include { make_db } from "./subworkflows/shotgun.nf"

workflow rnaseq {
    assemble()
}

workflow identify {
    search()
}

workflow combine_databases {
    make_db()
}
