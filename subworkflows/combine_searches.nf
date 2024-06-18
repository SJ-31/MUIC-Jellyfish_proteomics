include { COMBINE_ALL } from '../modules/combine_all.nf'
include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { COMBINE_PEP } from '../modules/combine_pep'
include { BLASTP } from '../modules/blastp'
include { SIGNALP } from '../modules/signalp'
include { DEEPLOC } from '../modules/deeploc'
include { SORT_BLAST } from '../modules/sort_blast'
include { ANNOTATE } from '../modules/annotate'
include { FINAL_METRICS } from '../modules/final_metrics'
include { MERGE_OPEN } from '../modules/merge_open'
include { CLUSTER_UNMATCHED } from '../modules/cluster_unmatched'
include { COVERAGE_CALC; COVERAGE_SPLIT; COVERAGE_MERGE } from '../modules/coverage'
include { INTERPROSCAN; SORT_INTERPRO } from '../modules/interpro'
include { EGGNOG; SORT_EGGNOG } from '../modules/eggnog'

def nullIfEmpty(ch, filename) {
    ch.ifEmpty({file("$params.storage/$filename")})
}

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
    blast_db

    main:
    SEARCH_INTERSECT(prot2intersect,
                     "$outdir/Combined", seq_header_mappings)
    COMBINE_PEP(psm2combinedPEP, true,
                    "$outdir/Combined")
    MERGE_OPEN(SEARCH_INTERSECT.out.unsorted, open_results,
                unmatched_pep_tsv,
               // 1. Merge combined database hits with results from open search
               // 2. Extract unmatched|denovo|transcriptome peptides as fasta file
        "$outdir/Combined")

    ranges = Channel.from(tuple(0, 35), tuple(35, 50),
                            tuple(50, 85), tuple(85, 100000))
    BLASTP(MERGE_OPEN.out.unknown_fasta.first(), blast_db,
        ranges, "$outdir/Unmatched/BLAST")
        // 1. Match de novo, transcriptome and unmatched peptides against proteins
        //     in database
        // * To optimize the searches, queries will be partitioned into four ranges,
        //  following the recommendations on the blast user guide
    BLASTP.out.collectFile(name: "combined_blast.csv", newLine: true)
        .set { combined_blast }

    SORT_BLAST(MERGE_OPEN.out.unknown_tsv,
                MERGE_OPEN.out.unmatched_pep,
                MERGE_OPEN.out.database_tsv,
                combined_blast, seq_header_mappings,
                MERGE_OPEN.out.unknown_fasta,
                "$outdir/Unmatched/BLAST")
    // 1. Merge blast results into identifications
    //      i.e. if a peptide A was matched to protein B, then A is added to the list of
    //      B's peptides in the main database search results
    // 2. Extract queries that weren't matched by blast

    EGGNOG(SORT_BLAST.out.unmatched,
            "$outdir/Unmatched/eggNOG")
    SORT_EGGNOG(EGGNOG.out.unmatched, // Extract peptides that weren't matched
                MERGE_OPEN.out.unmatched_pep, // by eggnog
                "$outdir/Unmatched/eggNOG")
    INTERPROSCAN(SORT_EGGNOG.out.fasta,
                    "$outdir/Unmatched/InterPro")
    SORT_INTERPRO(INTERPROSCAN.out, SORT_EGGNOG.out.unmatched,
                    unmatched_pep_tsv,
                    "$outdir/Unmatched/InterPro",
                    "$outdir/Unmatched/Remaining_unmatched")

    // BUG Temporary fix for eggnog failure
    eggnog_matched = nullIfEmpty(SORT_EGGNOG.out.matched, "EGGNOG_EMPTY")
    interpro_matched = nullIfEmpty(SORT_INTERPRO.out.matched, "INTERPRO_EMPTY")
    interpro_unmatched_fasta = nullIfEmpty(SORT_INTERPRO.out.unmatched_fasta,
        "INTERPRO_FASTA_EMPTY")

    ANNOTATE(SORT_BLAST.out.matched,
                "$outdir/Unmatched/Database-annotated")
    COMBINE_ALL(ANNOTATE.out.annotations, eggnog_matched,
                interpro_matched, directlfq, flashlfq, maxlfq,
                "$outdir", "$outdir/Logs")
    unmatched_ch = ANNOTATE.out.unannotated.mix(interpro_unmatched_fasta)
    // Key for "ProteinId" column:
    // P = protein from downloaded database
    // D = de novo peptide
    // T = transcriptome protein
    // U = peptide that was not matched to a protein

    COVERAGE_SPLIT(COMBINE_ALL.out.all)
    COVERAGE_CALC(COVERAGE_SPLIT.out.flatten(), "$outdir")
    COVERAGE_MERGE(COVERAGE_CALC.out.collect(), COMBINE_ALL.out.all, "$outdir")

    // Will not run if all proteins were matched
    CLUSTER_UNMATCHED(unmatched_ch.collect(), "$outdir")
    SIGNALP(CLUSTER_UNMATCHED.out.fasta, "$outdir/SignalP")
    DEEPLOC(SIGNALP.out.fasta, COVERAGE_MERGE.out.tsv,
        SEARCH_INTERSECT.out.unsorted,
        "$outdir/Deeploc")

}
