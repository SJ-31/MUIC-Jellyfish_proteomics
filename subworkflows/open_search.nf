include { MSFRAGGER as MSFRAGGER_OPEN } from "../modules/msfragger"
include { MSFRAGGER as MSFRAGGER_GLYCO } from "../modules/msfragger"
include { METAMORPHEUS as METAMORPHEUS_GPTMD } from "../modules/metamorpheus"
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
    MSFRAGGER_OPEN(mzML.collect(), "$params.config/open_params.params",
                   "GPTMD", "$outdir/MsFragger_open", "$outdir/Logs", dbWdecoys)
    MSFRAGGER_GLYCO(mzML.collect(), "$params.config/Nglyco_params.params",
                    "Glyco", "$outdir/MsFragger_glyco", "$outdir/Logs", dbWdecoys)
    METAMORPHEUS_GLYCO(mzML.collect(), "$outdir/Metamorpheus_glyco",
                       "$outdir/Logs",
                       "$params.config/metamorpheus_glyo_params.toml", "Glyco",
                       db)
    METAMORPHEUS_GPTMD(mzML.collect(), "$outdir/Metamorpheus_gptmd",
                       "$outdir/Logs",
                       "$params.config/metamorpheus_gptmd_params.toml", "GPTMD",
                       db)
}
