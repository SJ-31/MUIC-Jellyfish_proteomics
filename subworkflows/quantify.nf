include { FLASHLFQ } from '../modules/flashlfq'
include { MAP_SCANS } from '../modules/map_scans'
include { DIRECTLFQ; DIRECTLFQ_FORMAT } from '../modules/directlfq'
include { UNMATCHED_MSMS; UNMATCHED_PSMS } from '../modules/unmatched'


workflow 'quantify'{
    take:
    msms_mappings
    mzmls
    engine_percolator_output
    metamorpheus_AllPSMs
    maxquant_pin_file
    tide_target_search
    outdir

    main:
    engine_percolator_output.branch {
        comet: it =~ /comet/
        identipy: it =~ /identipy/
        msfragger: it =~ /msfragger/
        maxquant: it =~ /maxquant/
        tide: it =~ /tide/
        metamorpheus: it =~ /metamorpheus/
    }.set { per }
    MAP_SCANS(per.comet.mix(per.identipy, per.msfragger, maxquant_pin_file,
                            metamorpheus_AllPSMs, tide_target_search),
              msms_mappings,
              "$outdir/mapped_scans")
    FLASHLFQ(MAP_SCANS.out.collect(), mzmls, "$outdir")
    DIRECTLFQ_FORMAT(MAP_SCANS.out.collect(), msms_mappings, "$outdir")
    DIRECTLFQ(DIRECTLFQ_FORMAT.out, "$outdir")
    UNMATCHED_MSMS(per.comet, per.identipy, per.msfragger, maxquant_pin_file,
             metamorpheus_AllPSMs, tide_target_search,
                   msms_mappings, mzmls, "$outdir/Unmatched")
    UNMATCHED_PSMS(engine_percolator_output.collect(), "$outdir/Unmatched")
}
