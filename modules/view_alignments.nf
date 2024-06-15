process VIEW_ALIGNMENTS {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(identifications)
    path(alignments)
    val(outdir)
    //

    output:
    path("*.svg")
    //

    script:
    """
    view_alignments.py \
        -r $identifications \
        -c 0.8 \
        -a $alignments \
        -o .

    cp .command.out .log
    """
    //
}
