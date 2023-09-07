include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { ANNOTATE } from '../modules/annotate'
include { FINAL_METRICS } from '../modules/final_metrics'
include { MERGE_QUANT } from '../modules/merge_quantifications'
include { INTERPROSCAN } from '../modules/interpro'
include { UNMATCHED_PSMS } from '../modules/unmatched'

workflow 'combine_searches' {

    take:
    prot2intersect
    psm2combinedPEP
    directlfq
    percolator_psms
    outdir
    seq_header_mappings

    main:
    SEARCH_INTERSECT(prot2intersect,
                     "$outdir/Combined", seq_header_mappings)
    COMBINE_PEP(psm2combinedPEP, true,
                    "$outdir/Combined")
    MERGE_QUANT(directlfq, SEARCH_INTERSECT.out.sorted,
        "$outdir/Combined")
    UNMATCHED_PSMS(percolator_psms.collect(), "$outdir/Unmatched")

    ANNOTATE(MERGE_QUANT.out.database, "$outdir")
    // if ( params.denovo ) {
    //     BLASTP(MERGE_QUANT.out.unknown_fasta)
    //     INTERPROSCAN(, TODO: Give this the fasta file of peptides that don't
    //     have blast hits
    //     "$outdir")
    // } else {
    // FINAL_METRICS(MERGE_QUANT.out.database)
    // }
    //
    // Perform the emPAI calculations after de novo and transcriptome peptides have
    // been mapped back onto proteins


    emit:
    intersected_searches = SEARCH_INTERSECT.out.unsorted
    combinedPEP_psm = COMBINE_PEP.out
}
