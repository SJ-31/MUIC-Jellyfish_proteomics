Channel.fromPath(params.raws)
    .set { raw_ch }
Channel.fromPath(params.mzmls)
    .set { mzmls_ch }
Channel.fromPath(params.mgfs)
    .set { mgf_ch }




include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { MSGF } from '../modules/msgf'

workflow 'database_search' {
    MAXQUANT(raw_ch, params.mqpars, params.maxquant_dir)
    // MSFRAGGER(mzmls_ch, params.fragger_dir)
    // MSGF(mzmls_ch, params.msgf_dir)
}

workflow 'de_novo' {
    SMSNET(mgf_ch, params.smsnet_dir, params.smsnetmodel)
}
