process MAX_LFQ {
    publishDir "$outdir", mode: "copy"

    input:
    path(aqreformat)
    val(outdir)
    //

    output:
    path("max_lfq.tsv")
    //

    script:
    def check = file("${outdir}/max_lfq.tsv")
    if (check.exists()) {
        """
        cp ${outdir}/max_lfq.tsv .
        """
    } else {
        """
        Rscript ${params.bin}/R/iq_wrapper.r --input $aqreformat \
            --output max_lfq.tsv \
            --r_source ${params.bin}/R
        """
    }
    //
}
