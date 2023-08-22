/*
* initial checks and adapter filtering
*/
include { MULTIQC; TRIM; FASTQCREPORT; STATS } from '../modules/RNA-seq_initial_checks'
include { BBDUK; BLAST; SORTMERNA } from '../modules/RNA-seq_filtering'
include { GETFROMDIR; GRP_SEQKIT } from '../modules/Get_stats' addParams(ext: '.log', report_name: 'preprocessing_stats.txt' )

workflow clean {
    main:
    Channel
        .fromPath(params.rawfastq)
        .collect()
        .map { it -> ["0-initial_checks", it]}
        .set { raw_fastq_ch }

    // Get initial stats
    fastqc = FASTQCREPORT(raw_fastq_ch)
    Channel
        .fromFilePairs(params.rawpairs)
        .set { rawpairs_ch }

    filter_ch = Channel.fromPath(params.filters)
    overrepresented_ch = Channel.fromPath(params.queries)

    // Trim adapters
    fastp = TRIM(rawpairs_ch).flatten()
    fastp.branch{
        report: it =~ /html|json/
        fastqs: it =~ /fastq.gz/
        }
        .set { branched }
    branched.report.collect().map { it -> ["1-adapters_trimmed", it] }
        .set { trim_report } // For multiqc
    branched.fastqs.collect().map { it -> ["1-adapters_trimmed", it]}
        .set { grouped_fastqs } // For seqkit
    branched.fastqs.buffer( size: 2 ).map { it -> ["trimmed_only", it] }
        .set { trimmed }
    fastqc
        .mix(trim_report)
        .set { to_mqc } // Combine reports for multiqc
    MULTIQC(to_mqc)
    // Filter for rRNA
    BBDUK(trimmed, filter_ch).reads
        .set { bbduk_reads }
    SORTMERNA(trimmed, filter_ch.filter( ~/.*Chironex.*/ ))
        .set { srtmrna_results }
    STATS(raw_fastq_ch.mix(grouped_fastqs).mix(bbduk_reads)).collect()
        .set { summaries }
    bbduk_reads.flatten().filter(~/.*Chironex.*.fastq.gz/)
        .buffer(size: 2).map { it -> ["Chironex_filter", it ] }
        .set{ cx_filter }
    // Summarize stats
    GRP_SEQKIT(summaries)

    emit:
    trimmed
    cx_filter
}

/*
 * Assembly, assessment and translation
 */
include { SPADES; TRANSDECODE } from '../modules/RNA-seq_assembly'
include { GETSPADES } from '../modules/Get_stats.nf'

workflow assemble {
    main:
    clean()
    clean.out.trimmed.mix(clean.out.cx_filter)
        .set{ assemble_ch }
    GETSPADES(SPADES(assemble_ch, "$params.assembled"))
        .map { it -> ["rnaspades", it]}.set { spades_ch }
    TRANSDECODE(spades_ch, "$params.translated")
}
