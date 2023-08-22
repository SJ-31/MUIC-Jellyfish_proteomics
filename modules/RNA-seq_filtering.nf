params.k = "31"
params.cleaned = "../resources/RNA-seq/2-filtering"
params.logs = "$projectDir/logs"
params.normalized_seqs = "../resources/RNA-seq/3-normalize"
params.reports = "../reports/RNA-seq"
params.stats = "../reports/RNA-seq/2-filtering"
params.normalized = "../reports/RNA-seq/3-normalize"

process BBDUK {
    tag "Cleaning with $reference"
    memory { 2.GB * task.attempt }
    publishDir "$params.stats/${reference.baseName}", pattern: "*.txt", mode: 'copy'
    publishDir "$params.cleaned", pattern: "{*.fastq.gz,*-flagged*}"

    input:
    tuple val(id), path(reads) // You NEED the comma
    each path(reference)
    //
    output:
    tuple(val(reference.baseName), path("*.fastq.gz"), emit: reads)
    path("*-flagged*")
    path("*.txt")

    exec:
    def formatted = reads[0].baseName.replaceAll(/[12].*/, '') //.baseName converts paths into strs
    def forw = reads[0].baseName.replaceAll(/_T.*/, '')
    def rev = reads[1].baseName.replaceAll(/_T.*/, '')
    //
    script:
    """
    bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${reference.baseName}-${forw}_C.fastq.gz \
    out2=${reference.baseName}-${rev}_C.fastq.gz \
    ref=${reference} stats=${formatted}_${reference.baseName}.txt \
    outm=${formatted}-flagged.fastq.gz  k=$params.k -Xmx1G -Xms16M
    """
    // bbduk can run out of memory if you don't specify it
}

process SORTMERNA {
    tag "Cleaning with $reference"
    publishDir "$params.cleaned", pattern: "{*sortmerna*}"

    input:
    tuple val(id), path(reads)
    each path(reference)
    //
    output:
    path("*")

    exec:
    def formatted = reads[0].baseName.replaceAll(/[12].*/, '')
    //
    script:
    """
    sortmerna --reads ${reads[0]} --reads ${reads[1]} --ref ${reference} \
    --workdir sortmerna-${formatted}_${reference.baseName}
    """
}

process BLAST {
    publishDir "$params.cleaned", mode: 'copy'

    input:
    path(query)
    path(database)
    //
    output:
    path "*"
    //
    exec:
    def formatted = query.baseName.replaceAll(/_.*/, '') + '_BLAST-search.txt'

    script:
    """
    blastn -db $database -query $query -out $formatted
    """
    //
}
