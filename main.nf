include { search } from "./subworkflows/shotgun"
include { assemble } from "./subworkflows/rna_seq"
include { make_db } from "./subworkflows/make_db"
include { pre } from "./subworkflows/shotgun"
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow main.nf --manifest <manifest_file> --subworkflow <subworkflow>")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)


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

workflow {
    switch(params.subworkflow) {
        case "rnaseq":
            assemble();
            break;
        case "identify":
            println """
            Searching spectra...
                Prefix: $params.pref
                Manifest file: $params.manifest_file
                Database file: $params.db_spec
                Results path: $params.results
            """
            search(manifest.mzML, manifest.mgf, manifest.raw, db.normal);
            break;
        case "analyze":
            analysis(Channel.fromPath(params.first_final, params.second_final),
                "${params.results}/Analysis");
            break;
        case "preprocess":
            pre(manifest.raw, db.normal);
            break;
        case "combine_databases":
            println """
            Combining databases...
                Prefix: $params.pref
                Database file: $params.to_construct
                Manifest file: $params.manifest_file
                Results path: $params.results
            """
            if (params.to_construct == null) {
                println "Text file to databases is empty!"
                exit 0;
            }
            make_db(manifest.mzML, manifest.mgf)
            break;
    }
}

workflow hidden {

}
