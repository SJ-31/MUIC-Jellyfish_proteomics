process CONCAT_TSV {
    publishDir "$outdir", mode: 'copy'

    input:
    path(a)
    path(b)
    val(filename)
    val(outdir)
    //
    output:
    path("${filename}.tsv")

    //
    script:
    """
    awk '(NR == 1) || (FNR > 1)' ${a} ${b} > "${filename}.tsv"
    """
    //
}
