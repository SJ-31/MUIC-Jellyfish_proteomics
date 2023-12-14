include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { BLASTP } from '../modules/blastp'
include { SORT_BLAST } from '../modules/sort_blast'
include { ANNOTATE } from '../modules/annotate'
include { FINAL_METRICS } from '../modules/final_metrics'
include { MERGE_QUANT } from '../modules/merge_quantifications'
include { INTERPROSCAN } from '../modules/interpro'
include { EGGNOG; SORT_EGGNOG } from '../modules/eggnog'
include { UNMATCHED_PSMS } from '../modules/unmatched'

workflow 'combine_searches' {

    take:
    prot2intersect
    psm2combinedPEP
    directlfq
    percolator_psms
    outdir
    seq_header_mappings
    open_results

    main:
    SEARCH_INTERSECT(prot2intersect,
                     "$outdir/Combined", seq_header_mappings)
    COMBINE_PEP(psm2combinedPEP, true,
                    "$outdir/Combined")
    MERGE_QUANT(directlfq, SEARCH_INTERSECT.out.sorted,
        "$outdir/Combined")
    UNMATCHED_PSMS(percolator_psms.collect(), "$outdir/Unmatched")
    ANNOTATE(MERGE_QUANT.out.database_tsv, "$outdir")
    if ( params.denovo ) {
        // Syntax is <prefix> <do_one_hit?> <do_best_only?> <identity_threshold> <pep_threshold> <evalue_threshold>
        Channel.of("no_one_hits_no_degenerates 0 1 80 0.05 0.00001",
                   "one_hits_degenerates 1 0 80 0.05 0.00001",
                   "no_one_hits_degenerates 0 0 80 0.05 0.00001",
                   "one_hits_no_degenerates 1 1 80 0.05 0.00001")
            .set { blast_vars }

        BLASTP(MERGE_QUANT.out.unknown_fasta, params.blast_db,
               "$outdir/Unmatched/BLAST")
        SORT_BLAST(MERGE_QUANT.out.unknown_tsv, MERGE_QUANT.out.database_tsv,
                   BLASTP.out, seq_header_mappings, UNMATCHED_PSMS.out.fasta,
                   blast_vars, "$outdir/Unmatched") // Extract proteins that
                                                    // weren't matched by blast
        EGGNOG(SORT_BLAST.out.unmatched,
               "$outdir/Unmatched/eggNOG")
        SORT_EGGNOG(EGGNOG.out.unmatched, // Extract proteins that weren't matched
               "$outdir/Unmatched/eggNOG") // by eggnog
        INTERPROSCAN(SORT_EGGNOG.out,
                     "$outdir/Unmatched/InterPro")
        // FINAL_METRICS()
    } else {5
        FINAL_METRICS(MERGE_QUANT.out.database)
    }
    //
    // Perform the emPAI calculations after de novo and transcriptome peptides have
    // been mapped back onto proteins


    emit:
    intersected_searches = SEARCH_INTERSECT.out.unsorted
    combinedPEP_psm = COMBINE_PEP.out
}
