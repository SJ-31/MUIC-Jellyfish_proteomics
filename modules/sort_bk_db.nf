process MAKE_BK_DB {

    input:
    path(percolator_output)
    path(mapping_file)
    path(database)
    //

    output:
    tuple val(engine), path(sorted_database)
    //

    script:
    """

    """
    //
}
