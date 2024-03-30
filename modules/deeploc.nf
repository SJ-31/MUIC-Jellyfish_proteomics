process DEEPLOC {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/deeploc2"

    input:
    path(fasta)
    val(outdir)
    //

    output:
    path("deeploc_results.csv")
    //

    script:
    def check = file("deeploc_results.csv")
    if (check.exists()) {
    """
    mv -Z ${outdir}/deeploc_results.csv .
    """
    } else {
    """
    mkdir results
    deeploc2 -f $fasta \
            -o results
    cp results/* deeploc_results.csv
    """
    }
    //
}
