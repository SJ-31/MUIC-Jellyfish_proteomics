process PHILOSOPHER_FILTER {
    publishDir "$outdir", mode: "copy"

    input:
    path(pepxml)
    val(outdir)
    //

    output:
    //

    script:
    """
    philosopher filter \
        --pep 0.01

    """
    //
}
/* Options for FDR filtering
* --razor
* --cappedsequential
* --
*/
