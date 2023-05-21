// overwrite only gets rid of the files that have the same names...
params.rawpairs = "../resources/raw/RNA/*{1,2}.fastq"
params.queries = "../reports/*BLAST*"
params.rawfastq = "../resources/raw/RNA/*.fastq"
params.blastdb = "../resources/reference/unwanted_seqs/Chironex_BlASTDB/Chironex_DB"
params.filters = "../resources/reference/unwanted_seqs/*.fasta"
params.alignment = "../reports/RNA-seq/bowtie-read_alignment"
// params.reports = "$projectDir/report/RNA-seq" // If you don't specify $projectDir, then the directory will get created in nextflow's 'work' directory for that process
// Mix operator combines the items of several channels into one

Channel
    .fromPath(params.rawfastq, checkIfExists: true)
    .collect() // Collect groups items in the channel into a single list argument
    .map { it -> ["0-initial_checks", it]}
    .set { raw_fastq_ch }

Channel
    .fromFilePairs(params.rawpairs, checkIfExists: true)
    .set { rawpairs_ch }

filter_ch = Channel.fromPath(params.filters)
overrepresented_ch = Channel.fromPath(params.queries)

log.info """
    RNA-SEQ
    """.stripIndent()
println "Project : $workflow.projectDir"
println "Cmd line: $workflow.commandLine"
/*
* Initial checks and adapter filtering
*/
include { MULTIQC; TRIM; FASTQCREPORT; STATS; TODIR } from './modules/RNA-seq_initial_checks'
include { BBDUK; ORNA; BLAST; SORTMERNA } from './modules/RNA-seq_filtering'
include { GETFROMDIR; GRP_SEQKIT } from './modules/Get_stats' addParams(ext: '.log', report_name: 'preprocessing_stats.txt' )

workflow 'clean' {
    main:
    // Get initial stats
    fastqc = FASTQCREPORT(raw_fastq_ch)
    // Trim adapters
    fastp = TRIM(rawpairs_ch).flatten()
    fastp.branch{
        report: it =~ /html|json/
        fastqs: it =~ /fastq.gz/
        }
        .set { branched } // Split the channel into two
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
    GETFROMDIR(srtmrna_results)
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
include { BLOOM; SPADES; PLASS; TRANSDECODE; QALIGN; BALIGN } from './modules/RNA-seq_assembly'
include { GETBLOOM; GETSPADES } from './modules/Get_stats.nf'

workflow 'assemble' {
    main:
    clean() // We call the previous workflow so now we have access to the results from it
    clean.out.trimmed.mix(clean.out.cx_filter)
        .set{ assemble_ch }
    GETBLOOM(BLOOM(assemble_ch))
        .map { it -> ["rnabloom", it]}.set { bloom_ch }
    GETSPADES(SPADES(assemble_ch))
        .map { it -> ["rnaspades", it]}.set { spades_ch }
    TRANSDECODE(spades_ch.mix(bloom_ch))
    spades_ch.mix(bloom_ch)
        .map { it -> [ (it[1].baseName =~ /(.*)_[A-Z]/)[0][1], it[0], it[1] ]
        }
        .join(assemble_ch.mix(assemble_ch))
        .set { reads_assemblies_ch } // Need to call on the channel twice so that won't miss out on the other species
    PLASS(assemble_ch)
    // spades_ch.mix(bloom_ch)
    QALIGN(reads_assemblies_ch) // The alignment metrics aren't working
    reads_assemblies_ch
    BALIGN(reads_assemblies_ch, params.alignment)

    // emit:
}

workflow {
    assemble()

}


