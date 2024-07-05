process VIEW_ALIGNMENTS {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(identifications)
    path(alignments)
    path(peptide_map)
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
        -o . \
        -m "engine_alignment" \
        -p $peptide_map

    cp .command.out .log
    """
    //
}
