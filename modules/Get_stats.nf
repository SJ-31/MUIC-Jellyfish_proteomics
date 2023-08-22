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
