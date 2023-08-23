process DEISOTOPE {
    publishDir "$outdir", mode: "copy"

    input:
    path(mzMLs)
    path(outdir)
    //

    output:
    path("*")
    //

    script:
    """
    deisotope.sh $file_list .
    """
    //
}
