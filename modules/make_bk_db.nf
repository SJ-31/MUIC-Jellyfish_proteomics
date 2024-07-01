process MAKE_BK_DB {
    publishDir "$outdir", mode: "copy"

    input:
    tuple val(engine), path(percolator_valid), path(percolator_decoys)
    path(database)
    val(outdir)
    //

    output:
    tuple val(engine), path(output)
    //

    script:
    output = "${engine}_bk_database.fasta"
    if (file("${outdir}/${output}").exists()) {
        """
        cp "${outdir}/${output}" .
        """
    } else {
        """
        make_bk_db.py $percolator_valid \
            $percolator_decoys \
            $database \
            "${engine}_bk_database.fasta"
        wc -l "${engine}_bk_database.fasta"
        grep ">" "${engine}_bk_database.fasta" | wc -l
        """
    }
    //
}
