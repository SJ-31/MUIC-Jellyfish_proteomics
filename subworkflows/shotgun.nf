// CONVENTION: the final protein/PSMs output file will be prefixed with the engine that produced them
include { MAXQUANT } from '../modules/maxquant'
include { MSFRAGGER } from '../modules/msfragger'
include { SMSNET } from '../modules/smsnet'
include { COMET } from '../modules/comet'
include { IDENTIPY } from '../modules/identipy'
include { METAMORPHEUS } from '../modules/metamorpheus'
include { TIDE } from '../modules/tide'
include { TIDE_COMBINED_PEP } from '../modules/tide'
include { CASANOVO } from '../modules/casanovo'
include { PERCOLATOR } from '../modules/percolator'
include { MS2RESCORE } from '../modules/ms2rescore'
include { COMBINED_DATABASE } from '../modules/combined_database'
include { EXTRACT_CASANOVO } from '../modules/casanovo'
include { DEISOTOPE } from '../modules/deisotope'
include { PEPNET } from '../modules/pepnet'
include { EXTRACT_PEPNET } from '../modules/pepnet'
include { bk_decoys } from './bk_decoys.nf'
include { combine_searches as combine_searches_FIRST } from './combine_searches.nf'
include { combine_searches as combine_searches_SECOND } from './combine_searches.nf'

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
    if ( params.denovo ) {
    SMSNET(manifest.mgf.collect(), "$params.results/SMSNET")
    EXTRACT_CASANOVO(
            CASANOVO(manifest.mzML.collect(),"$params.results/Casanovo"),
            "$params.results/Casanovo")
    PEPNET(manifest.mgf, "$params.results/PepNet")
    EXTRACT_PEPNET(PEPNET.out.collect(), "$params.results/PepNet")
    EXTRACT_CASANOVO.out.mix(SMSNET.out, EXTRACT_PEPNET.out).flatten()
    .branch {
        combined: it ~/.*combined.*/
        decoys: it ~/.*decoys.*/
        normal: it ~/.*normal.*/
        }.set { denovo }
    } else {
        denovo = Channel.empty()
    }
    COMBINED_DATABASE(database_listing, denovo.combined,
    denovo.normal, denovo.decoys,
        "$projectDir/data/reference/protein_databases/combined")
}

dbWdecoys = params.databaseWdecoy
db = params.database
workflow 'search' {
    empty = Channel.empty()
    // MaxQuant seems to only work with .raw files
    // MAXQUANT(manifest.raw, "$params.results/MaxQuant", db)
    COMET(manifest.mzXML.collect(), "$params.results/First_pass/Comet",
    dbWdecoys)
    MSFRAGGER(manifest.mzML.collect(), "$params.config/MSFragger_params.params",
    "$params.results/First_pass/MsFragger", dbWdecoys)
    IDENTIPY(manifest.mzML.collect(), "$params.results/First_pass/Identipy", dbWdecoys)
    METAMORPHEUS(manifest.mzML.collect(), "$params.results/First_pass/Metamorpheus", db)
    // MSGF(manifest.mzML, "$params.results/First_pass/msgf", db) No longer used,
    //  no way to integrate with percolator for now
    TIDE(manifest.mzXML.collect(), "$params.results/First_pass/Tide", "$params.results/First_pass/Percolator",
    db)
    TIDE_COMBINED_PEP(TIDE.out.percolator, "$params.results/First_pass/Percolator")
    TIDE.out.percolator.flatten().filter( ~/.*\.target\.proteins\.txt/ )
        .set { tide_percolator }
    // MAXQUANT.out.ms2rescore.flatten().filter( ~/.*txt/ ).collect()
    //     .map { it -> ["maxquant", it] }
    //     .set { maxqms2rescore }
    // MS2RESCORE(maxqms2rescore, "$params.results/First_pass/MaxQuant",
    // manifest.mgf.collect())
    empty.mix(
        METAMORPHEUS.out.percolator,
        //MS2RESCORE.out,
        COMET.out.percolator,
        MSFRAGGER.out.percolator,
        IDENTIPY.out.percolator
    ).set { to_percolator }


    PERCOLATOR(to_percolator, "$params.results/First_pass/Percolator", dbWdecoys)

    combine_searches_FIRST(PERCOLATOR.out.prot2intersect.mix(tide_percolator),
                           PERCOLATOR.out.psm2combinedPEP
                                .mix(TIDE_COMBINED_PEP.out.psm2combinedPEP),
                           PERCOLATOR.out.prot2combinedPEP
                                .mix(TIDE_COMBINED_PEP.out.prot2combinedPEP),
                        "$params.results/First_pass/Combined")
    // Second pass with Bern and Kil decoy database
    bk_decoys(PERCOLATOR.out.prot, dbWdecoys, manifest.mzXML, manifest.mzML)
    // PERCOLATOR.out.filter( !()
}
