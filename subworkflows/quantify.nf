include { UNMATCHED_MSMS } from '../modules/unmatched'
include { MAP_SCANS } from '../modules/map_scans'
include { DIRECTLFQ_FORMAT } from '../modules/directlfq'


workflow 'quantify'{
    take:
    msms_mappings
    mzmls
    percolator_psms
    percolator_proteins
    metamorpheus_AllPSMs
    tide_target_search
    psm2combinedPEP
    outdir

    main:
    percolator_psms.filter{ !(it =~ /.*metamorpheus.*|.*tide.*/) }
        .mix(metamorpheus_AllPSMs, tide_target_search)
        .map { it = [ it.baseName.replaceFirst(/_.*/, ""), it ] }
        .set { psm }
    percolator_proteins
        .map { it = [ it.baseName.replaceFirst(/_.*/, ""), it ] }
        .set { prot }

    MAP_SCANS(psm.join(prot),
              msms_mappings,
              "$outdir/Mapped_scans")
    DIRECTLFQ_FORMAT(MAP_SCANS.out.collect(), msms_mappings, "$outdir")
    UNMATCHED_MSMS(MAP_SCANS.out.collect(), psm2combinedPEP.collect(), mzmls,
                "$outdir/Unmatched")

    emit:
    unmatched_msms = UNMATCHED_MSMS.out.mzmls
    directlfq_input = DIRECTLFQ_FORMAT.out

}
