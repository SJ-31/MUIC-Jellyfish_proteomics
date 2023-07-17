process GRP {
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
