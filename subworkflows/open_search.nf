
include { MSFRAGGER as MSFRAGGER_OPEN } from "../modules/msfragger"
include { MSFRAGGER as MSFRAGGER_GLYCO } from "../modules/msfragger"
include { METAMORPHEUS as METAMORPHEUS_GET_GPTMD } from "../modules/metamorpheus"
include { PERCOLATOR } from "../modules/percolator"
include { MAP_SCANS } from '../modules/map_scans'
include { MS_MAPPING } from '../modules/ms_mapping'
include { METAMORPHEUS as METAMORPHEUS_SEARCH_GPTMD } from "../modules/metamorpheus"
include { METAMORPHEUS as METAMORPHEUS_GLYCO } from "../modules/metamorpheus"
include { SORT_OPEN } from "../modules/sort_open_searches"


workflow 'open_search' {
    // Special database searches on the unmatched ms/ms spectra, or those that did not pass a certain threshold
    //
    take:
    mzML
    dbPlusDecoys
    db
    outdir
    seq_header_mapping

    main:
    MS_MAPPING(mzML.collect(), "$outdir")
    MSFRAGGER_OPEN(mzML.collect(), "$params.config_dir/open_fragger.params",
                   "GPTMD", "$outdir/MsFragger_open", "$outdir/Logs", dbPlusDecoys)
    MSFRAGGER_GLYCO(mzML.collect(), "$params.config_dir/Nglyco_fragger.params",
                    "Glyco", "$outdir/MsFragger_glyco", "$outdir/Logs",
                    dbPlusDecoys)
    METAMORPHEUS_GET_GPTMD(mzML.collect(), "$outdir/Metamorpheus_gptmd",
                       "$outdir/Logs", "",
                       "$params.config_dir/metamorpheus_gptmd_params.toml", db)
    METAMORPHEUS_GET_GPTMD.out.all.flatten().filter( ~/.*\.xml/ )
        .set { xml_ch }
    METAMORPHEUS_SEARCH_GPTMD(mzML.collect(),
                              "$outdir/Metamorpheus_gptmd_search",
                              "$outdir/Logs", "GTPMD",
                              "$params.config_dir/metamorpheus_params.toml", xml_ch)
    METAMORPHEUS_SEARCH_GPTMD.out.percolator.mix(
        MSFRAGGER_OPEN.out.percolator,
        MSFRAGGER_GLYCO.out.percolator,
    ).set { to_percolator }
    PERCOLATOR(to_percolator, "$outdir/Percolator",
               "$outdir/Logs", dbPlusDecoys.first())
    SORT_OPEN(PERCOLATOR.out.prot2intersect.collect(),
              seq_header_mapping,
              "$outdir")

    PERCOLATOR.out.psms.filter({ !(it =~ /.*metamorpheus.*/ )})
                            .mix(METAMORPHEUS_SEARCH_GPTMD.out.psms)
                            .map { it = [ it.baseName.replaceFirst(/_.*/, ""), it ] }
                            .set { psm }

   PERCOLATOR.out.prot2intersect
        .map { it = [ it.baseName.replaceFirst(/_.*/, ""), it ] }
        .set { prot }

    MAP_SCANS(psm.join(prot),
              Channel.fromPath("$params.config/NO_FILE"),
              MS_MAPPING.out,
              "$outdir/Mapped_scans")

    emit:
    open_results = SORT_OPEN.out
    prot2intersect = PERCOLATOR.out.prot2intersect
}
