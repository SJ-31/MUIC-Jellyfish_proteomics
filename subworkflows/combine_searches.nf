include { COMBINE_ALL } from '../modules/combine_all.nf'
include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { BLASTP } from '../modules/blastp'
include { SORT_BLAST } from '../modules/sort_blast'
include { ANNOTATE } from '../modules/annotate'
include { FINAL_METRICS } from '../modules/final_metrics'
include { MERGE_OPEN } from '../modules/merge_open'
include { COVERAGE_CALC; COVERAGE_SPLIT; COVERAGE_MERGE } from '../modules/coverage'
include { INTERPROSCAN; SORT_INTERPRO } from '../modules/interpro'
include { EGGNOG; SORT_EGGNOG } from '../modules/eggnog'

workflow 'combine_searches' {

    take:
    prot2intersect
    psm2combinedPEP
    flashlfq
    maxlfq
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
        ranges = Channel.from(tuple(0, 35), tuple(35, 50),
                              tuple(50, 85), tuple(85, 100000))
        BLASTP(MERGE_OPEN.out.unknown_fasta.first(),
               params.blast_db, ranges,
               "$outdir/Unmatched/BLAST") // Let blast try to annotate
        // de novo, transcriptome and unmatched peptides
        // To optimize the searches, queries will be partitioned into four ranges,
        // following the recommendations on the blast user guide
        BLASTP.out.collectFile(name: "combined_blast.csv", newLine: true)
            .set { combined_blast }

        SORT_BLAST(MERGE_OPEN.out.unknown_tsv,
                   MERGE_OPEN.out.unmatched_pep,
                   MERGE_OPEN.out.database_tsv,
                   combined_blast, seq_header_mappings,
                   MERGE_OPEN.out.unknown_fasta,
                   "$outdir/Unmatched/BLAST") // Extract peptides that
                                                    // weren't matched by blast
        EGGNOG(SORT_BLAST.out.unmatched,
               "$outdir/Unmatched/eggNOG")
        SORT_EGGNOG(EGGNOG.out.unmatched, // Extract peptides that weren't matched
                    MERGE_OPEN.out.unmatched_pep, // by eggnog
                    "$outdir/Unmatched/eggNOG")
        INTERPROSCAN(SORT_EGGNOG.out.fasta,
                     "$outdir/Unmatched/InterPro")
        SORT_INTERPRO(INTERPROSCAN.out, SORT_EGGNOG.out.unmatched,
                      "$outdir/Unmatched/InterPro",
                      "$outdir/Unmatched/Remaining_unmatched")
        ANNOTATE(SORT_BLAST.out.matched,
                 "$outdir/Unmatched/Database-annotated")
        COMBINE_ALL(ANNOTATE.out.annotations, SORT_EGGNOG.out.matched,
                    SORT_INTERPRO.out.matched,
                    directlfq, flashlfq, maxlfq, "$outdir")
        // Key for "ProteinId" column:
        // P = protein from downloaded database
        // D = de novo peptide
        // T = transcriptome protein
        // U = peptide that was not matched to a protein
    } else {
        Channel.fromPath("$params.config_dir/NO_FILE")
            .map{ it = ["no_file", it, it] }
            .set { no_file }
        ANNOTATE(MERGE_OPEN.out.database_tsv, "$outdir")
        COMBINE_ALL(no_file, no_file, ANNOTATE.out.annotations,
                    directlfq, flashlfq, maxlfq, "$outdir")
    }

    COVERAGE_SPLIT(COMBINE_ALL.out.all)
    COVERAGE_CALC(COVERAGE_SPLIT.out.flatten(), "$outdir")
    COVERAGE_MERGE(COVERAGE_CALC.out.collect(), COMBINE_ALL.out.all,
                   "$outdir")

    emit:
    all_combined = Channel.empty()
    all_combined = COMBINE_ALL.out.all
    intersected_searches = SEARCH_INTERSECT.out.unsorted
    combinedPEP_psm = COMBINE_PEP.out
}
