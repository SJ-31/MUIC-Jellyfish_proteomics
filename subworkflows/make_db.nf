include { CASANOVO } from '../modules/casanovo'
include { EXTRACT_CASANOVO } from '../modules/casanovo'
include { COMBINED_DATABASE } from '../modules/combined_database'
include { PEPNET } from '../modules/pepnet'
include { EXTRACT_PEPNET } from '../modules/pepnet'

Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.indexed_mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mgf: it =~ /.mgf/
        mzXML: it =~ /.mzXML/
        raw: it =~ /.raw/
    }.set { manifest }

workflow 'make_db' {

    main:
    Channel.fromPath(params.to_construct)
        .splitText()
        .set { database_listing }
    if ( params.denovo ) {
    // SMSNET(mgf.collect(), "$params.results/SMSNET") // TODO: Fixthis
    CASANOVO(manifest.mzML,"$params.results/Casanovo")
        EXTRACT_CASANOVO(CASANOVO.out.peps.collect(), "$params.results/Casanovo")
    PEPNET(manifest.mgf, "$params.results/PepNet")
    EXTRACT_PEPNET(PEPNET.out.peps.collect(), "$params.results/PepNet")

    EXTRACT_CASANOVO.out.mix(EXTRACT_PEPNET.out)
     .set { denovo }
    } else {
        denovo = Channel.empty()
    }
    COMBINED_DATABASE(database_listing.collect(), denovo.collect(),
                      "$params.results/databases")

    emit:
    COMBINED_DATABASE.out.listing
}
