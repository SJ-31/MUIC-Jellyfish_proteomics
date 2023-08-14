process DEISOTOPE {
    publishDir "$outdir", mode: "copy"

    input:
    path(file_list)
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
