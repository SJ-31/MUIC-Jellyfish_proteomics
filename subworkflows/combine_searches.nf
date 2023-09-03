include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { ANNOTATE } from '../modules/annotate'
include { MERGE_QUANT } from '../modules/merge_quantifications'
include { INTERPROSCAN } from '../modules/interpro'

workflow 'combine_searches' {

    take:
    prot2intersect
    psm2combinedPEP
    directlfq
    outdir
    header_mappings
    seq_mappings

    main:
    SEARCH_INTERSECT(prot2intersect,
                     "$outdir/Combined", header_mappings, seq_mappings)
    COMBINE_PEP(psm2combinedPEP, true,
                    "$outdir/Combined")
    MERGE_QUANT(directlfq, SEARCH_INTERSECT.out.sorted, seq_mappings,
        "$outdir/Combined")
    ANNOTATE(MERGE_QUANT.out, seq_mappings, "$outdir")
    if ( params.denovo ) {
        INTERPROSCAN(ANNOTATE.out.denovo.mix(ANNOTATE.out.transcriptome),
                    "$outdir")
    }




    emit:
    intersected_searches = SEARCH_INTERSECT.out.unsorted
    combinedPEP_psm = COMBINE_PEP.out
}
