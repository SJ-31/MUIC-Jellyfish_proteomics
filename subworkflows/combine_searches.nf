include { COMBINE_ALL } from '../modules/combine_all.nf'
include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { BLASTP } from '../modules/blastp'
include { SORT_BLAST } from '../modules/sort_blast'
include { ANNOTATE } from '../modules/annotate'
include { FINAL_METRICS } from '../modules/final_metrics'
include { MERGE_OPEN } from '../modules/merge_open'
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
    MERGE_OPEN(SEARCH_INTERSECT.out.sorted, open_results,
                unmatched_pep_tsv, // Merge combined database hits
                // with quantification, and obtain unmatched peptides
        "$outdir/Combined")
    if (params.denovo) {
        // Syntax is <prefix> <do_one_hit?> <do_best_only?> <identity_threshold> <pep_threshold> <evalue_threshold>
        BLASTP(MERGE_OPEN.out.unknown_fasta, params.blast_db,
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
        SORT_BLAST(MERGE_OPEN.out.unknown_tsv,
                   MERGE_OPEN.out.unmatched_pep,
                   MERGE_OPEN.out.database_tsv,
                   BLASTP.out, seq_header_mappings,
                   MERGE_OPEN.out.unknown_fasta,
                   blast_vars, "$outdir/Unmatched/BLAST") // Extract peptides that
                                                    // weren't matched by blast
        EGGNOG(SORT_BLAST.out.unmatched,
               "$outdir/Unmatched/eggNOG")
        SORT_EGGNOG(EGGNOG.out.unmatched, // Extract peptides that weren't matched
                    MERGE_OPEN.out.unmatched_pep.first(), // by eggnog
                    "$outdir/Unmatched/eggNOG")
        INTERPROSCAN(SORT_EGGNOG.out.fasta,
                     "$outdir/Unmatched/InterPro")
        SORT_INTERPRO(INTERPROSCAN.out, SORT_EGGNOG.out.unmatched,
                      "$outdir/Unmatched/InterPro",
                      "$outdir/Unmatched/Remaining_unmatched")
        ANNOTATE(SORT_BLAST.out.matched,
                 "$outdir/Unmatched/Database-annotated")

        ANNOTATE.out.annotations
            .join(SORT_EGGNOG.out.matched)
            .join(SORT_INTERPRO.out.matched)
            .flatten()
            .branch {
                eggnog: it[1] ==~ /.*eggnog.*/
                interpro: it[1] ==~ /.*interpro.*/
                download: it[1] ==~ /.*download.*/
            }. set { combined }

        combined.eggnog.view()
        combined.interpro.view()
        combined.download.view()

        // COMBINE_ALL(combine.eggnog, combine.interpro, combine.download,
        // directlfq, flashlfq, "$outdir")
    } else {
        Channel.fromPath("$params.config/NO_FILE")
            .map{ it = ["no_file", it, it] }
            .set { no_file }
        ANNOTATE(MERGE_OPEN.out.database_tsv, "$outdir")
        COMBINE_ALL(no_file, no_file, ANNOTATE.out.annotations,
                    directlfq, flashlfq, "$outdir")
    }

    emit:
    intersected_searches = SEARCH_INTERSECT.out.unsorted
    combinedPEP_psm = COMBINE_PEP.out
}
