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

    main:
    SEARCH_INTERSECT(prot2intersect.collect(),
                     "$outdir/Combined")
    COMBINE_PEP_PSM(psm2combinedPEP.collect(), true,
                    "$outdir/Combined")
    COMBINE_PEP_PROT(prot2combinedPEP.collect(), false,
                    "$outdir/Combined")
    ANNOTATE(SEARCH_INTERSECT.out, "$params.results/First_pass")

    emit:
    intersected_searches = SEARCH_INTERSECT.out
    combinedPEP_psm = COMBINE_PEP_PSM.out
    combinedPEP_prot = COMBINE_PEP_PROT.out
}
