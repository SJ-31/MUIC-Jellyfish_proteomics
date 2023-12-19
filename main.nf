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

Channel.fromPath(params.db_spec).splitText() { it.replaceAll("\n", "") }
    .branch {
        normal: it ==~ /.*all_normal.fasta/
        plusdecoys: it ==~ /.*decoysWnormal.fasta/
        seq_mapping : it ==~ /.*decoysWnormal_mapping.tsv/
        header_mapping : it ==~ /.*header_mappings.tsv/
        downloaded : it ==~ /.*downloaded.fasta/
    }.set { db }

workflow rnaseq {
    assemble()
}

workflow identify {
    search(manifest.mzML, manifest.mgf, manifest.raw, db.normal)
}

workflow preprocess {
    pre(manifest.mzML, manifest.raw, db.normal)
}

workflow combine_databases {
    make_db(manifest.mzML, manifest.mgf)
}

workflow hidden {

}
