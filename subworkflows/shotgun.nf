Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mzXML: it =~ /.mzXML/
        mgf: it =~ /.mgf/
        raw: it =~ /.raw/
    }.set { manifest }

include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { MSGF } from '../modules/msgf'
include { COMET } from '../modules/comet'
include { IDENTIPY } from '../modules/identipy'
include { METAMORPHEUS } from '../modules/metamorpheus'
include { CASANOVO } from '../modules/casanovo'
include { PERCOLATOR } from '../modules/percolator'

workflow 'search' {
    // MaxQuant seems to only work with .raw files
    // MAXQUANT(manifest.raw, params.mqpars, "$params.results/MaxQuant")
    // .set { maxq }
    // COMET(manifest.mzXML.collect(), "$params.results/Comet")
    // .set { comet }
    // MSFRAGGER(manifest.mzML.collect(), params.fragger_closedPars, "$params.results/MsFragger")
    // .set { fragger }
    IDENTIPY(manifest.mzML.collect(), "$params.results/Identipy")
        .set { ipy }
    // METAMORPHEUS(manifest.mzML.collect(), "$params.results/Metamorpheus")
    //     .set{ mmorph }
    // MSGF(manifest.mzML.collect(), "$params.results/msgf")
    // SMSNET(manifest.mgf.collect(), "$params.results/SMSNET")
    // CASANOVO(manifest.mzML.collect(), "$params.results/Casanovo")
    // PERCOLATOR(
    //         .mix(
    //             mmorph.percolator,
    //             maxq.percolator,
    //             comet.percolator,
    //             fragger.percolator,
    //             ipy.percolator
    // )

}
