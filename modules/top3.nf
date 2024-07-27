process TOP3 {
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_results)
    path(aqreformat)
    val(outdir)
    //

    output:
    path("top3.tsv")
    //

    script:
    check = "top3.tsv"
    if (file("${outdir}/${check}").exists()) {
        """
        cp ${outdir}/top3.tsv .
        """
    } else {
    """
    helpers.py -t top3 \
        -i $combined_results \
        -d $aqreformat \
        -o top3.tsv
    """
    }
    //
}
