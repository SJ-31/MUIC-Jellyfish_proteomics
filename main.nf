include { search } from "./subworkflows/shotgun"
include { analyze } from "./subworkflows/analysis"
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


Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.indexed_mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mgf: it =~ /.mgf/
        raw: it =~ /.raw/
    }.set { manifest }

def getDatabases() {
    db = Channel.fromPath(params.db_spec).splitText() { it.replaceAll("\n", "") }
        .branch {
            normal: it ==~ /.*all_normal.fasta/
            plusdecoys: it ==~ /.*decoysWnormal.fasta/
            seq_header_mapping : it ==~ /.*mappings.tsv/
            downloaded : it ==~ /.*downloaded.fasta/
        }
    return db
}

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
            db = getDatabases()
            search(manifest.mzML, manifest.mgf, db.normal, db.plusdecoys,
                   db.seq_header_mapping);
            break;
        case "analyze":
            first_final = "${params.results}/1-First_pass/${params.pref}_all_wcoverage.tsv"
            alignments = "${params.results}/1-First_pass/aligned_peptides.tsv"
            peptide_map = "${params.results}/1-First_pass/percolator_peptide_map.tsv"
            analyze(Channel.fromPath(first_final),
                 Channel.fromPath(alignments),
                 Channel.fromPath(peptide_map),
                "${params.results}/Analysis");
            break;
        case "preprocess":
            pre(manifest.raw, manifest.mzML, getDatabases().normal);
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
