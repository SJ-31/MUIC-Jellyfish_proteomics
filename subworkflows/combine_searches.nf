include { COMBINE_ALL } from '../modules/combine_all.nf'
include { SEARCH_INTERSECT } from '../modules/search_intersect'
include { BLASTP } from '../modules/blastp'
include { SIGNALP } from '../modules/signalp'
include { DEEPLOC } from '../modules/deeploc'
include { SORT_BLAST } from '../modules/sort_blast'
include { ANNOTATE } from '../modules/annotate'
include { WRITE_QUANT; LFQ_MERGE } from '../modules/write_quant'
include { TOP3 } from '../modules/top3'
include { CONCAT_TSV } from '../modules/utils'
include { DIRECTLFQ } from '../modules/directlfq'
include { MAX_LFQ } from '../modules/maxlfq'
include { FASTS } from '../modules/fasts'
include { MERGE_OPEN } from '../modules/merge_open'
include { CLUSTER_UNMATCHED } from '../modules/cluster_unmatched'
include { COMBINE_PERCOLATOR } from '../modules/combine_percolator'
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
    outdir
    seq_header_mappings
    open_results
    prot2intersect_open_search
    directlfq_input
    blast_db

    main:
    SEARCH_INTERSECT(prot2intersect,
                     "$outdir/Combined", seq_header_mappings)
    MERGE_OPEN(SEARCH_INTERSECT.out.unsorted, open_results,
               // 1. Merge combined database hits with results from open search
               // 2. Extract denovo|transcriptome peptides as fasta file
        "$outdir/Combined")


    ranges = Channel.from(tuple(0, 35), tuple(35, 50),
                            tuple(50, 85), tuple(85, 100000))
    BLASTP(MERGE_OPEN.out.unknown_fasta.first(), blast_db,
        ranges, "$outdir/Unmatched/BLAST")
        // 1. Match de novo and transcriptome against proteins
        //     in database
        // * To optimize the searches, queries will be partitioned into four ranges,
        //  following the recommendations on the blast user guide
    BLASTP.out.collectFile(name: "combined_blast.csv", newLine: true)
        .set { combined_blast }
    queries = MERGE_OPEN.out.unknown_fasta_sep.flatten()

    FASTS(queries, params.fasts_db,
        "$outdir/Unmatched/BLAST")
    FASTS.out.collectFile(name: "combined_fasts.csv", newLine: true)
        .set { combined_fasts }

    SORT_BLAST(MERGE_OPEN.out.unknown,
                MERGE_OPEN.out.database_tsv,
                combined_blast, combined_fasts,
                seq_header_mappings,
                "$outdir/Unmatched/BLAST")
    // 1. Merge blast results into identifications
    //      i.e. if a peptide A was matched to protein B, then A is added to the list of
    //      B's peptides in the main database search results
    // 2. Extract queries that weren't matched by blast

    // Optional annotation
    if (params.annotate) {
        EGGNOG(SORT_BLAST.out.unmatched,
                "$outdir/Unmatched/eggNOG")
        SORT_EGGNOG(EGGNOG.out.unmatched, // Extract peptides that weren't matched by eggnog
                    "$outdir/Unmatched/eggNOG")
        INTERPROSCAN(SORT_EGGNOG.out.fasta,
                        "$outdir/Unmatched/InterPro")
        SORT_INTERPRO(INTERPROSCAN.out, SORT_EGGNOG.out.unmatched,
                        "$outdir/Unmatched/InterPro",
                        "$outdir/Unmatched/Remaining_unmatched")

        // BUG Temporary fix for eggnog failure
        eggnog_matched = nullIfEmpty(SORT_EGGNOG.out.matched, "EGGNOG_EMPTY")
        interpro_matched = nullIfEmpty(SORT_INTERPRO.out.matched, "INTERPRO_EMPTY")
        interpro_unmatched_fasta = nullIfEmpty(SORT_INTERPRO.out.unmatched_fasta,
            "INTERPRO_FASTA_EMPTY")


        ANNOTATE(SORT_BLAST.out.matched,
                    "$outdir/Unmatched/Database-annotated")

        ANNOTATE.out.annotations
            .concat(eggnog_matched, interpro_matched).collect()
            .set { combined_ch }
    } else {
        SORT_BLAST.out.matched.concat(
            nullIfEmpty(Channel.empty(), "EGGNOG_EMPTY"),
            nullIfEmpty(Channel.empty(), "INTERPRO_EMPTY")
            ).collect().set { combined_ch }
    }

    COMBINE_ALL(combined_ch, "$outdir", "$outdir/Logs")
    COMBINE_PERCOLATOR(prot2intersect, prot2intersect_open_search,
                        COMBINE_ALL.out.all,
                        seq_header_mappings,
                        "$outdir", "$outdir/Logs")

    WRITE_QUANT(COMBINE_ALL.out.all, directlfq_input, "$outdir/Quantify")
    TOP3(COMBINE_ALL.out.all, directlfq_input, "$outdir/Quantify")
    DIRECTLFQ(WRITE_QUANT.out, "$outdir/Quantify", "$outdir/Logs")
    MAX_LFQ(WRITE_QUANT.out, "$outdir/Quantify")

    lfq = DIRECTLFQ.out.quant.concat(TOP3.out, MAX_LFQ.out).collect()
    LFQ_MERGE(lfq, seq_header_mappings, "$outdir")

    COVERAGE_SPLIT(COMBINE_ALL.out.all)
    COVERAGE_CALC(COVERAGE_SPLIT.out.flatten(),
        COMBINE_PERCOLATOR.out.seq_map.first(),
        COMBINE_PERCOLATOR.out.peptide_map.first(),
        "$outdir/.coverage")
    COVERAGE_MERGE(COVERAGE_CALC.out.collect(), COMBINE_ALL.out.all, "$outdir")


    // Will not run if all proteins were matched
    if (params.annotate) {
        unmatched_ch = ANNOTATE.out.unannotated.mix(interpro_unmatched_fasta)
        CLUSTER_UNMATCHED(unmatched_ch.collect(), "$outdir")
        SIGNALP(CLUSTER_UNMATCHED.out.fasta, "$outdir/SignalP")
        DEEPLOC(SIGNALP.out.fasta, COVERAGE_MERGE.out.tsv,
            SEARCH_INTERSECT.out.unsorted,
            "$outdir/Deeploc", "$outdir")
    }
}
