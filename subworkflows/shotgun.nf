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
include { COMBINED_DATABASE } from '../modules/combined_database'
include { EXTRACT_CASANOVO } from '../modules/extract_casanovo'
include { DEISOTOPE } from '../modules/deisotope'
include { bk_decoys } from './bk_decoys.nf'

workflow 'preprocess' {
    DEISOTOPE(params.manifest_file)
}

Channel.fromPath(params.manifest_file)
    .splitCsv(header: true, sep: "\t")
    .map { it -> [ it.Prefix, it.Raw, it.mzML, it.mzXML, it.mgf ] }
    .flatten().branch {
        mzML: it =~ /.mzML/
        mzXML: it =~ /.mzXML/
        mgf: it =~ /.mgf/
        raw: it =~ /.raw/
    }.set { manifest }

// Channel.fromPath(params.databases)
//     .splitText()
//     .set { database_listing }



workflow 'make_db' {
    // if ( params.denovo ) {
    // SMSNET(manifest.mgf.collect(), "$params.results/SMSNET")
      //   .set { smsnet_ch }
    // EXTRACT_CASANOVO(CASANOVO(manifest.mzML.collect(),"$params.results/Casanovo"),
    // "$params.results/Casanovo")
    //     .set { casanovo_ch }
//     casanovo_ch.mix(smsnet_ch).flatten()
//     .branch {
    //     combined: it ~/*combined*/
    //     decoys: it ~/*decoys*/
    //     normal: it ~/*normal*/
    //     }.set { denovo }
    // } else {
    //     denovo = Channel.empty()
    // }
    // COMBINED_DATABASE(database_listing, denovo.combined,
    // denovo.normal, denovo.decoys,
        // "$projectDir/data/reference/protein_databases/combined")
}

dbWdecoys = params.databaseWdecoy
db = params.database
workflow 'search' {
    empty = Channel.empty()
    // MaxQuant seems to only work with .raw files
    // MAXQUANT(manifest.raw, "$params.results/MaxQuant", db)
    COMET(manifest.mzXML.collect(), "$params.results/Comet",
    dbWdecoys)
    MSFRAGGER(manifest.mzML.collect(), "$params.config/MSFragger_params.params",
    "$params.results/MsFragger", dbWdecoys)
    IDENTIPY(manifest.mzML.collect(), "$params.results/Identipy", dbWdecoys)
    METAMORPHEUS(manifest.mzML.collect(), "$params.results/Metamorpheus", db)
    MSGF(manifest.mzML, "$params.results/msgf", db)
    TIDE(manifest.mzXML.collect(), "$params.results/Tide", "$params.results/Percolator",
    db)
    TIDE.out.percolator.flatten().filter( ~/.*\.target\.proteins\.txt/ )
        .set { tide_percolator }
    // MAXQUANT.out.ms2rescore.flatten().filter( ~/.*txt/ ).collect()
    //     .map { it -> ["maxquant", it] }
    //     .set { maxqms2rescore }
    // MS2RESCORE(maxqms2rescore, "$params.results/MaxQuant",
    // manifest.mgf.collect())
    empty.mix(
        METAMORPHEUS.out.percolator,
        //MS2RESCORE.out,
        COMET.out.percolator,
        MSFRAGGER.out.percolator,
        IDENTIPY.out.percolator
    ).set { to_percolator }

    PERCOLATOR(to_percolator, "$params.results/Percolator", dbWdecoys)
    SEARCH_INTERSECT(PERCOLATOR.out.prot2intersect.mix(tide_percolator).collect(),
                     "$params.results/Combined")
    COMBINE_PEP_PSM(PERCOLATOR.out.psm2combinedPEP.collect(), true,
                    "$params.results/Combined")
    COMBINE_PEP_PROT(PERCOLATOR.out.prot2combinedPEP.collect(), false,
                    "$params.results/Combined")

    // Second pass with Bern and Kil decoy database
    bk_decoys(PERCOLATOR.out.prot, dbWdecoys, manifest.mzXML, manifest.mzML)
    // bk_decoys.out.prot2intersect.view()

}
