process METAMORPHEUS {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/metamorpheus"

    input:
    path(raws)
    val(outdir)
    //

    output:
    //

    script:
    """
    metamorpheus -s $raws \
        -o . \
        -t $params.morpheusPars \
        -d $params.databaseWdecoy \
    """
    // Need to add extraction for percolator
}
