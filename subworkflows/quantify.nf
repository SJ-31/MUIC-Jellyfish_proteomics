include { FLASHLFQ } from '../modules/flashlfq'
include { MAP_SCANS } from '../modules/map_scans'
include { DIRECTLFQ; DIRECTLFQ_FORMAT } from '../modules/directlfq'
include { FILTER_MSMS; UNMATCHED_PSMS } from '../modules/unmatched'


workflow 'quantify'{
    take:
    msms_mappings
    mzmls
    engine_percolator_output
    metamorpheus_AllPSMs
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
    MAP_SCANS(per.comet.mix(per.identipy, per.msfragger, per.maxquant,
                            metamorpheus_AllPSMs, tide_target_search),
              msms_mappings,
              "../$outdir/Mapped_scans")
    FLASHLFQ(MAP_SCANS.out.collect(), mzmls.collect(), "$outdir")
    DIRECTLFQ_FORMAT(MAP_SCANS.out.collect(), msms_mappings, "$outdir")
    DIRECTLFQ(DIRECTLFQ_FORMAT.out, "$outdir")

    MAP_SCANS.out
        .map { it -> [ it.baseName.replaceAll("_.*", ""), it ] }
        .set { mapped_scans }
    engine_percolator_output
        .map { it -> [ it.baseName.replaceAll("_.*", ""), it ] }
        .set { perc_psms }
    mapped_scans.join(perc_psms)
        .set { scans_psms }
    FILTER_MSMS(scans_psms, msms_mappings, mzmls, "../$outdir/Unmatched")
    UNMATCHED_PSMS(engine_percolator_output.collect(), "../$outdir/Unmatched")

    emit:
    directlfq = DIRECTLFQ.out


}
