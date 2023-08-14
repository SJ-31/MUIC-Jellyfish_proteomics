process MAKE_BK_DB {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(engine), path(percolator_valid), path(percolator_decoys)
    path(mapping_file)
    path(database)
    val(outdir)
    //

    output:
    tuple val(engine), path("${engine}_bk_database.fasta")
    //

    script:
    """
    make_bk_db.py $percolator_valid \
        $percolator_decoys \
        $database \
        $mapping_file \
        "${engine}_bk_database.fasta"
    """
    //
}
