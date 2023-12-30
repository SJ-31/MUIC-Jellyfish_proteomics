include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { BLASTP } from '../modules/blastp'
include { SORT_BLAST } from '../modules/sort_blast'
include { ANNOTATE } from '../modules/annotate'
include { FINAL_METRICS } from '../modules/final_metrics'
include { MERGE_QUANT } from '../modules/merge_quantifications'
include { INTERPROSCAN; SORT_INTERPRO } from '../modules/interpro'
include { EGGNOG; SORT_EGGNOG } from '../modules/eggnog'

workflow 'combine_searches' {

    take:
    prot2intersect
    psm2combinedPEP
    flashlfq
    directlfq
    unmatched_pep_tsv
    outdir
    seq_header_mappings
    open_results

    main:
    SEARCH_INTERSECT(prot2intersect,
                     "$outdir/Combined", seq_header_mappings)
    COMBINE_PEP(psm2combinedPEP, true,
                    "$outdir/Combined")
    MERGE_QUANT(directlfq, flashlfq,
                SEARCH_INTERSECT.out.sorted, open_results,
                unmatched_pep_tsv, // Merge combined database hits
                // with quantification, and obtain unmatched peptides
        "$outdir/Combined")
    if (params.denovo) {
        // Syntax is <prefix> <do_one_hit?> <do_best_only?> <identity_threshold> <pep_threshold> <evalue_threshold>
        BLASTP(MERGE_QUANT.out.unknown_fasta, params.blast_db,
               "$outdir/Unmatched/BLAST") // Let blast try to annotate
        // de novo, transcriptome and unmatched peptides
        Channel.of("nO_nD 0 1 80 0.05 0.00001",
                   "O_D 1 0 80 0.05 0.00001",
                   "nO_D 0 0 80 0.05 0.00001",
                   "O_nD 1 1 80 0.05 0.00001")
            .set { blast_vars }
        // nO_nD = No one hits, no degenerates:
        //      targets matched by one query only not allowed
        //      peptides that match multiple targets not allowed
        // O_D = One hits, degenerates
        // nO_D = No one, hits, degenerates
        // O_nD
        SORT_BLAST(MERGE_QUANT.out.unknown_tsv,
                   MERGE_QUANT.out.unmatched_pep,
                   MERGE_QUANT.out.database_tsv,
                   BLASTP.out, seq_header_mappings,
                   MERGE_QUANT.out.unknown_fasta,
                   blast_vars, "$outdir/Unmatched/BLAST") // Extract peptides that
                                                    // weren't matched by blast
        EGGNOG(SORT_BLAST.out.unmatched,
               "$outdir/Unmatched/eggNOG")
        SORT_EGGNOG(EGGNOG.out.unmatched, // Extract peptides that weren't matched
                    MERGE_QUANT.out.unmatched_pep.first(), // by eggnog
                    "$outdir/Unmatched/eggNOG")
        INTERPROSCAN(SORT_EGGNOG.out.fasta,
                     "$outdir/Unmatched/InterPro")
        SORT_INTERPRO(INTERPROSCAN.out, SORT_EGGNOG.out.unmatched,
                      "$outdir/Unmatched/InterPro",
                      "$outdir/Unmatched/Remaining_unmatched")
        ANNOTATE(SORT_BLAST.out.matched,
                 "$outdir/Unmatched/Database-annotated")
        // TODO need to merge this with the kept blast results
        // combine_all(annotate.out, )
        // FINAL_METRICS()
    } else {
        ANNOTATE(MERGE_QUANT.out.database_tsv, "$outdir")
        FINAL_METRICS(MERGE_QUANT.out.database)
    }
    //
    // Perform the emPAI calculations after de novo and transcriptome peptides have
    // been mapped back onto proteins


    emit:
    intersected_searches = SEARCH_INTERSECT.out.unsorted
    combinedPEP_psm = COMBINE_PEP.out
}
