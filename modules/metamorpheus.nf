process METAMORPHEUS {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/metamorpheus"

    input:
    path(raws)
    val(outdir)
    //

    output:

    path("${params.pref}-comet_percolator")
    //

    script:
    """
    metamorpheus -s $raws \
        -o . \
        -t $params.morpheusPars \
        -d $params.databaseWdecoy \
    mv Task1SearchTask/{All*,results.txt} .
    percolator_wrapper.sh \
        -p $params.pref \
        -i AllPSMS_FormattedForPercolator.tab \
        -f $params.databaseWdecoy \
        -e metamorpheus
    """
    // Need to add extraction for percolator
}
