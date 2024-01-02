include { FLASHLFQ } from '../modules/flashlfq'
include { UNMATCHED_PSMS } from '../modules/unmatched'
include { MAP_SCANS } from '../modules/map_scans'
include { DIRECTLFQ; DIRECTLFQ_FORMAT } from '../modules/directlfq'
include { FILTER_MSMS } from '../modules/unmatched'


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

    UNMATCHED_PSMS(percolator_psms.collect(), "$outdir/Unmatched")
    MAP_SCANS(psm.join(prot),
              UNMATCHED_PSMS.out.tsv,
              msms_mappings,
              "$outdir/Mapped_scans")
    FLASHLFQ(MAP_SCANS.out.collect(), mzmls.collect(), "$outdir", "$outdir/Logs")
    DIRECTLFQ_FORMAT(MAP_SCANS.out.collect(), msms_mappings, "$outdir")
    DIRECTLFQ(DIRECTLFQ_FORMAT.out, "$outdir", "$outdir/Logs")
    FILTER_MSMS(MAP_SCANS.out.collect(), psm2combinedPEP.collect(), mzmls,
                "$outdir/Unmatched")

    emit:
    directlfq = DIRECTLFQ.out.quant
    flashlfq = FLASHLFQ.out.prot
    unmatched_msms = FILTER_MSMS.out
    unmatched_pep_tsv = UNMATCHED_PSMS.out.tsv

}
