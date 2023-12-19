include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { BLASTP } from '../modules/blastp'
include { SORT_BLAST } from '../modules/sort_blast'
include { ANNOTATE } from '../modules/annotate'
include { FINAL_METRICS } from '../modules/final_metrics'
include { MERGE_QUANT } from '../modules/merge_quantifications'
include { INTERPROSCAN; SORT_INTERPRO } from '../modules/interpro'
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
    UNMATCHED_PSMS(percolator_psms.collect(), "$outdir/Unmatched")
    MERGE_QUANT(directlfq, SEARCH_INTERSECT.out.sorted, open_results,
                UNMATCHED_PSMS.out.fasta, // Merge combined database hits
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
                   UNMATCHED_PSMS.out.tsv,
                   MERGE_QUANT.out.database_tsv,
                   BLASTP.out, seq_header_mappings,
                   MERGE_QUANT.out.unknown_fasta,
                   blast_vars, "$outdir/Unmatched") // Extract peptides that
                                                    // weren't matched by blast
        EGGNOG(SORT_BLAST.out.unmatched,
               "$outdir/Unmatched/eggNOG")
        SORT_EGGNOG(EGGNOG.out.unmatched, // Extract peptides that weren't matched
                    UNMATCHED_PSMS.out.tsv, // by eggnog
                    "$outdir/Unmatched/eggNOG")
        INTERPROSCAN(SORT_EGGNOG.out.fasta,
                     "$outdir/Unmatched/InterPro")
        SORT_INTERPRO(INTERPROSCAN.out, SORT_EGGNOG.out.unmatched,
                      "$outdir/Unmatched/InterPro")
        // ANNOTATE(SORT_BLAST.out.matched, "$outdir") TODO need to merge this with the kept blast results
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
