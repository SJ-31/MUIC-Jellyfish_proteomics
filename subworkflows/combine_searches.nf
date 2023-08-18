include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP as COMBINE_PEP_PSM } from '../modules/combine_pep'
include { COMBINE_PEP as COMBINE_PEP_PROT } from '../modules/combine_pep'
include { ANNOTATE } from '../modules/annotate'

workflow 'combine_searches' {

    take:
    prot2intersect
    psm2combinedPEP
    prot2combinedPEP
    outdir
    header_mappings

    main:
    SEARCH_INTERSECT(prot2intersect,
                     "$outdir/Combined", header_mappings)
    COMBINE_PEP_PSM(psm2combinedPEP, true,
                    "$outdir/Combined")
    COMBINE_PEP_PROT(prot2combinedPEP, false,
                    "$outdir/Combined")
    ANNOTATE(SEARCH_INTERSECT.out.unsorted, "$outdir")

    emit:
    intersected_searches = SEARCH_INTERSECT.out.unsorted
    combinedPEP_psm = COMBINE_PEP_PSM.out
    combinedPEP_prot = COMBINE_PEP_PROT.out
}
