include { CASANOVO; EXTRACT_CASANOVO } from '../modules/casanovo'
include { COMBINED_DATABASE } from '../modules/combined_database'
include { PEPNET; EXTRACT_PEPNET } from '../modules/pepnet'
include { SMSNET; EXTRACT_SMSNET } from '../modules/smsnet'


workflow 'make_db' {

    take:
    mzML
    mgf

    main:
    Channel.fromPath(params.to_construct)
        .splitText()
        .set { database_listing }
    if ( params.denovo ) {
        CASANOVO(mzML,"$params.results/Denovo/Casanovo")
            EXTRACT_CASANOVO(CASANOVO.out.peps.collect(),
                             "$params.results/Denovo/Casanovo")
        PEPNET(mgf, "$params.results/Denovo/PepNet")
        EXTRACT_PEPNET(PEPNET.out.peps.collect(), "$params.results/Denovo/PepNet")
        SMSNET(mgf, "$params.results/Denovo/SMSNet")
        EXTRACT_SMSNET(SMSNET.out.collect(), "$params.results/Denovo/SMSNet")

        EXTRACT_CASANOVO.out.mix(EXTRACT_PEPNET.out, EXTRACT_SMSNET.out)
        .set { denovo }
    } else {
        denovo = Channel.fromPath("$params.config_dir/empty.fasta")
    }
    COMBINED_DATABASE(database_listing.collect(), denovo.collect(),
                      "$params.results/Databases", params.db_loc)
}
