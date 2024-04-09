process RAWPARSE {
    publishDir "$outdir", mode: "copy"

    input:
    path(rawfiles)
    val(outdir)
    //

    output:
    path("*.mzML")
    //

    script:
    def check = file("${outdir}/${rawfiles[0].baseName}.mzML")
    if (check.exists()) {
        """
        cp ${outdir}/*.mzML .
        """
    } else {
        """
        mkdir raw; mv $rawfiles raw
        mono $params.thermo_parser -d raw \
            -f 2
        mv raw/*.mzML .
        """
    }
    // initially had a flag to disable peak picking, but this prevented the mzML
    // files from being used in the downstream functions
}
