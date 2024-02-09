include { FLASHLFQ } from '../modules/flashlfq'
include { UNMATCHED_PSMS; UNMATCHED_MSMS } from '../modules/unmatched'
include { MAP_SCANS } from '../modules/map_scans'
include { DIRECTLFQ; DIRECTLFQ_FORMAT } from '../modules/directlfq'
include { WRITE_QUANT } from '../modules/write_quant'


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
    UNMATCHED_MSMS(MAP_SCANS.out.collect(), psm2combinedPEP.collect(), mzmls,
                "$outdir/Unmatched")
    WRITE_QUANT(DIRECTLFQ.out.quant, FLASHLFQ.out.prot, "$outdir")

    emit:
    directlfq = WRITE_QUANT.out.dlfq
    flashlfq = WRITE_QUANT.out.flfq
    unmatched_msms = UNMATCHED_MSMS.out
    unmatched_pep_tsv = UNMATCHED_PSMS.out.tsv

}
