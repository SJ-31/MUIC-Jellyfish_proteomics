include { search } from "./subworkflows/shotgun"
include { assemble } from "./subworkflows/rna_seq"
include { make_db } from "./subworkflows/make_db"

workflow rnaseq {
    assemble()
}

workflow identify {
    search(Channel.from("placeholder"))
}

workflow combine_databases {
    make_db()
}

workflow full_pipeline {
    make_db() | search()
}
