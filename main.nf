include { search } from "./subworkflows/shotgun"
include { assemble } from "./subworkflows/rna_seq"
include { make_db } from "./subworkflows/make_db"
include { pre } from "./subworkflows/shotgun"

Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.indexed_mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mgf: it =~ /.mgf/
        raw: it =~ /.raw/
    }.set { manifest }

workflow rnaseq {
    assemble()
}

workflow identify {
    search(manifest.mzML, manifest.mgf, manifest.raw,
           Channel.from("placeholder"))
}

workflow preprocess {
    pre(manifest.mzML)
}

workflow combine_databases {
    make_db(manifest.mzML, manifest.mgf)
}

workflow full_pipeline {
    make_db() | search()
}
