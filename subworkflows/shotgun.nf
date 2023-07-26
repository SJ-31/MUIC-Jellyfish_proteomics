Channel.fromPath(params.raws)
    .set { raw_ch }
Channel.fromPath(params.mzmls)
    .set { mzmls_ch }
Channel.fromPath(params.mzXMLs)
    .set { mzxml_ch }
Channel.fromPath(params.mgfs)
    .set { mgf_ch }


include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { MSGF } from '../modules/msgf'
include { COMET } from '../modules/comet'
include { IDENTIPY } from '../modules/identipy'

workflow 'search' {
    // MAXQUANT(raw_ch, params.mqpars, "$params.results/MaxQuant")
    COMET(mzxml_ch, "$params.results/Comet")
    // MSFRAGGER(mzmls_ch, "$params.results/MsFragger")
    IDENTIPY(mzmls_ch, params.database,"$params.results/Identipy")
    // MSGF(mzmls_ch, "$params.results/msgf")
    // SMSNET(mgf_ch, "$params.results/SMSNET")
}
