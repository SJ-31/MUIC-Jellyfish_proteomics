include { MSFRAGGER as MSFRAGGER_OPEN } from "../modules/msfragger"
include { MSFRAGGER as MSFRAGGER_GLYCO } from "../modules/msfragger"
include { METAMORPHEUS as METAMORPHEUS_GET_GPTMD } from "../modules/metamorpheus"
include { PERCOLATOR } from "../modules/percolator"
include { METAMORPHEUS as METAMORPHEUS_SEARCH_GPTMD } from "../modules/metamorpheus"
include { METAMORPHEUS as METAMORPHEUS_GLYCO } from "../modules/metamorpheus"


workflow 'open_search' {
    // Special database searches on the unmatched ms/ms spectra, or those that did not pass a certain threshold
    //
    take:
    mzML
    dbWdecoys
    db
    outdir

    main:
    MSFRAGGER_OPEN(mzML.collect(), "$params.config/open_fragger.params",
                   "GPTMD", "$outdir/MsFragger_open", "$outdir/Logs", dbWdecoys)
    MSFRAGGER_GLYCO(mzML.collect(), "$params.config/Nglyco_fragger.params",
                    "Glyco", "$outdir/MsFragger_glyco", "$outdir/Logs",
                    dbWdecoys)
    METAMORPHEUS_GET_GPTMD(mzML.collect(), "$outdir/Metamorpheus_gptmd",
                       "$outdir/Logs", "",
                       "$params.config/metamorpheus_gptmd_params.toml", db)
    METAMORPHEUS_GET_GPTMD.out.all.flatten().filter( ~/.*\.xml/ )
        .set { xml_ch }
    METAMORPHEUS_SEARCH_GPTMD(mzML.collect(),
                              "$outdir/Metamorpheus_gptmd",
                              "$outdir/Logs", "GTPMD",
                              "$params.config/metamorpheus_params.toml", xml_ch)
    METAMORPHEUS_SEARCH_GPTMD.out.percolator.mix(
        MSFRAGGER_OPEN.out.percolator,
        MSFRAGGER_GLYCO.out.percolator,
    ).set { to_percolator }
}
