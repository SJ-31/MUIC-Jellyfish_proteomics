params.rawfastq = "../resources/raw/RNA/*.fastq"
params.stats = "../reports/RNA-seq/2-filtering"
params.reports = "../reports/RNA-seq" // If you don't specify $projectDir, then the directory will get created in nextflow's 'work' directory for that process for an output process
params.rnadata = "../resources/RNA-seq"
// Mix operator combines the items of several channels into one

process FASTQCREPORT {
    // tag "Report in $params.reports"

    publishDir "$params.reports/$loc", mode: 'copy'  // Only files matching the declaration in output: get published
    input:
    tuple val(loc), path(fastq)
    //
    output:
    tuple val("0-initial_checks"), path("fastqc")
    //
    script:// Also, you need to surround parameters with curly braces to reference them in shell commands
    """
    if test -f "fastqc"
    then
        fastqc $fastq -o fastqc -q
    else
        mkdir fastqc
        fastqc $fastq -o fastqc -q
    fi
    """
}

process MULTIQC {
    publishDir "$params.reports/$loc", mode: 'copy'

    input:
    tuple val(loc), path("*")
    //
    output:
    path "*"
    //
    script:
    """
    multiqc .
    """
    //
}

process TRIM {
    // tag "Report in $params.reports/1-adapters_trimmed"
    publishDir "$params.rnadata/1-adapters_trimmed", pattern: "*T.fastq.gz"
    publishDir "$params.reports/1-adapters_trimmed", pattern: "{*.html,*.json}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    //
    output:
    tuple val('trim'), path("*")

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} --detect_adapter_for_pe -c \
    -o ${reads[0].baseName}_T.fastq.gz -O ${reads[1].baseName}_T.fastq.gz -R ${sample_id}_trimming \
    -l 50 -z 4 -h ${sample_id}_trim.html -j ${sample_id}fastp.json
    """
    //
}

process TODIR {
    input:
    tuple val(id), path('*')
    //
    output:
    tuple val(id), path('*')
    //
    exec:

    script:
    """
    mkdir $id
    find . -type f | xargs -I{} mv {} $id
    """
    //
}

process STATS {
    publishDir "$params.reports/$loc", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(loc), path(files)
    //
    output:
    path '*.txt'
    //
    exec:
    formatted = loc.replaceAll(/.*-/, '')
    if ( loc =~ /C/ ) {
        id = files[0].baseName.replaceAll(/.*-/, '').replaceAll(/[12\.].*/, '')
        loc = "2-filtering/${formatted}"
        formatted = "${formatted}-${id}"
    }

    script:
    """
    seqkit stats $files -a -b -T > ${formatted}_seqkit-stats.txt
    """
    //
}
