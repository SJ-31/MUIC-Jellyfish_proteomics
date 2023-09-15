include { FLASHLFQ } from '../modules/flashlfq'
include { MAP_SCANS } from '../modules/map_scans'
include { DIRECTLFQ; DIRECTLFQ_FORMAT } from '../modules/directlfq'
include { FILTER_MSMS } from '../modules/unmatched'


workflow 'quantify'{
    take:
    msms_mappings
    mzmls
    engine_percolator_output
    metamorpheus_AllPSMs
    tide_target_search
    psm2combinedPEP
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
    MAP_SCANS(per.comet.mix(per.identipy, per.msfragger, per.maxquant,
                            metamorpheus_AllPSMs, tide_target_search),
              msms_mappings,
              "$outdir/Mapped_scans")
    FLASHLFQ(MAP_SCANS.out.collect(), mzmls.collect(), "$outdir", "$outdir/Logs")
    DIRECTLFQ_FORMAT(MAP_SCANS.out.collect(), msms_mappings, "$outdir")
    DIRECTLFQ(DIRECTLFQ_FORMAT.out, "$outdir", "$outdir/Logs")
    FILTER_MSMS(MAP_SCANS.out.collect(), psm2combinedPEP.collect(), mzmls,
                "$outdir/Unmatched")

    emit:
    directlfq = DIRECTLFQ.out.quant
    unmatched_msms = FILTER_MSMS.out



}
