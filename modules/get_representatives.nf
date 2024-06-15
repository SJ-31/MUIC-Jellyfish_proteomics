process GET_REPRESENTATIVES {
    publishDir "$outdir", mode: "copy"
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("${params.pref}_all_representatives.tsv")
    path("*.log")
    //

    script:
    """
    clustering.py \
        -m $params.mmseqs \
        -i $combined_tsv \
        -o "${params.pref}_all_representatives.tsv"

    cp .command.out .log
    """
    //
}
