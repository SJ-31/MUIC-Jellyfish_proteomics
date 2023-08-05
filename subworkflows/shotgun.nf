Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mzXML: it =~ /.mzXML/
        mgf: it =~ /.mgf/
        raw: it =~ /.raw/
    }.set { manifest }

// CONVENTION: the final protein/PSMs output file will be prefixed with the engine that produced them
include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { MSGF } from '../modules/msgf'
include { COMET } from '../modules/comet'
include { IDENTIPY } from '../modules/identipy'
include { METAMORPHEUS } from '../modules/metamorpheus'
include { TIDE } from '../modules/tide'
include { CASANOVO } from '../modules/casanovo'
include { PERCOLATOR } from '../modules/percolator'
include { MS2RESCORE } from '../modules/ms2rescore'
include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP as COMBINE_PEP_PSM } from '../modules/combine_pep'
include { COMBINE_PEP as COMBINE_PEP_PROT } from '../modules/combine_pep'

workflow 'denovo' {
    // SMSNET(manifest.mgf.collect(), "$params.results/SMSNET")
    // CASANOVO(manifest.mzML.collect(), "$params.results/Casanovo")
}

workflow 'make_db' {

}

dbWdecoys = $params.databaseWdecoy
db = $params.database
workflow 'search' {
    // MaxQuant seems to only work with .raw files
    MAXQUANT(manifest.raw, "$params.results/MaxQuant", db)
        .set { maxq }
    // COMET(manifest.mzXML.collect(), "$params.results/Comet",
    // dbWdecoys)
    // .set { comet }
    // MSFRAGGER(manifest.mzML.collect(), "$params.config/MSFragger_params.params",
    // "$params.results/MsFragger", dbWdecoys)
    // .set { fragger }
    // IDENTIPY(manifest.mzML.collect(), "$params.results/Identipy", dbWdecoys)
    //     .set { ipy }
    // METAMORPHEUS(manifest.mzML.collect(), "$params.results/Metamorpheus", db)
    //     .set { mmorph }
    // MSGF(manifest.mzML.collect(), "$params.results/msgf", db)
    //     .set { msgf }
    // TIDE(manifest.mzXML.collect(), "$params.results/Tide", "$params.results/Percolator",
    // db)
    //     .set { tide }
    MS2RESCORE(maxq.ms2rescore.collect(), "$params.results/MaxQuant",
    manifest.mgf.collect())
        .set { maxq_percolator }
    // PERCOLATOR(
    //         .mix(
    //             mmorph.percolator,
    //             maxq.percolator,
    //             comet.percolator,
    //             fragger.percolator,
    //             ipy.percolator,
    //             maxq_percolator,
    //             ),
    // "$params.results/Percolator")
    //     .set { percolator }
    // SEARCH_INTERSECT(percolator.prot2intersect.collect(), "$params.results/Combined")
    // COMBINE_PEP_PSM(percolator.psm2combinedPEP.collect(), true,
    //                 "$params.results/Combined")
    // COMBINE_PEP_PSM(percolator.prot2combinedPEP.collect(), false,
    //                 "$params.results/Combined")
}
