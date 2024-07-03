process MAKETREE {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(alignment)
    val(tool)
    val(outdir)
    //

    output:
    path("*.log")
    //

    script:
    if (tool == "iqtree2") {
        """
        iqtree2 -s ${alignment}

        """
    }

    //
}
