params.assembled = "../resources/RNA-seq/4-assembly"
params.stats = "../reports/RNA-seq/2-filtering"
params.want = "log"
params.ext = ""

process GETFROMDIR {
    publishDir params.stats, mode: 'copy', pattern: '*params.ext'

    input:
    path(folder)
    //
    output:
    path("*${params.want}*")
    //
    exec:
    formatted = "${folder.baseName}${params.ext}"

    script:
    """
    find . -name "*${params.want}*" | xargs -I{} cp {} ./$formatted
    """
    //
}

params.report_name = ''
process GRP_SEQKIT {
    publishDir params.stats, mode: 'copy'

    input:
    path(files)

    output:
    path(params.report_name)

    script: // You need the single quote separator for shell scripts
    """
    awk 'FNR==1 && NR!=1{next;}{print}' ${files} > $params.report_name
    """
}

process GETBLOOM {
    publishDir "$params.assembled/rnabloom"

    input:
    path(assembly)
    //
    output:
    path(formatted)
    //
    exec:
    formatted = "${assembly.baseName}.fasta"

    shell:
    '''
    cat $(find -L . -name *.fa) > !{formatted}
    '''
    //
}

process GETSPADES {
    publishDir "$params.assembled/rnaspades"

    input:
    path(assembly)
    //
    output:
    path(formatted)
    //
    exec:
    formatted = "${assembly.baseName}.fasta"

    shell:
    '''
    cp $(find -L . -maxdepth 2 -name transcripts.fasta) ./!{formatted}
    '''
    //
}

process A_STATS {
    input:

    //
    output:

    //
    script:
    """

    """
    //
}
