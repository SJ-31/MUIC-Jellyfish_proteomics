Channel.fromPath(params.mq_peptides).collect()
    .set { mq_peptides_ch }
Channel.fromPath(params.mzmls)
    .set { mzmls_ch }
Channel.fromPath(params.mgfs)
    .set { mgf_ch }

par_ch = Channel.fromPath(params.mqpars)


include { MAXQUANT } from './modules/maxquant'
include { MSFRAGGER } from './modules/msfragger'
include { SMSNET } from './modules/smsnet'

workflow 'database_search' {
    // MAXQUANT(par_ch)
    MSFRAGGER(mzmls_ch, params.fragger_dir)
}

workflow 'de_novo' {
    SMSNET(mgf_ch, params.smsnet_dir, params.smsnetmodel)
}

workflow 'collect_results' {
    Channel.fromPath(params.mq_peptides)
        .collectFile(name: "MaxQuant_AllPeptides.txt", keepHeader: true,
                     skip: 1, storeDir: "../reports/Proteomics")
        .set { mq_peptides_ch }
}

workflow {
    database_search()
    collect_results()
    de_novo()
}
