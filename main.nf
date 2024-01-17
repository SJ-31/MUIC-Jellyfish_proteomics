include { search } from "./subworkflows/shotgun"
include { assemble } from "./subworkflows/rna_seq"
include { make_db } from "./subworkflows/make_db"
include { pre } from "./subworkflows/shotgun"

/* LAUNCH FLAGS */
/* --entry <rnaseq|identify|analyze|preprocess|combine_databases> */
/*      rnaseq: runs transcriptome assembly pipeline */
/*      preprocess: preprocess and calibrate mass spectrometry files */
/*          this expects raw
/*      combine_databases: combine fasta files to create decoy database, */
/*          optionally run de novo search engines /*
/* --manifest_file <FILE> */
/*      tsv file containing paths to raw, mgf and indexed mzML files generated */
/*          from the preprocessing stage */
/*          two headers are "indexed_mzML", "mgf" and "raw"
/* --to_construct <FILE> */
/*      text file containing all the protein sources for */
/*          the database search engines, in fasta format */
/* --denovo */
/*      whether or not to run de novo search engines on the mass
/*          spec files */



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

workflow analyze {
    analysis(Channel.fromPath(params.first_final, params.second_final),
             "${params.results}/Analysis")
}

workflow preprocess {
    pre(manifest.mzML, manifest.raw, db.normal)
}

workflow combine_databases {
    make_db(manifest.mzML, manifest.mgf)
}

workflow hidden {

}
