process INTERPROSCAN {
    publishDir "$outdir", mode: "copy"

    input:
    path(unknown_fasta)
    val(outdir)
    //

    output:
    path("${unknown_fasta.baseName}-SCAN.tsv")
    //

    script:
    """
    interproscan.sh -i $unknown_fasta \
        -f tsv \
        -goterms \
        -pa \
        -b ${unknown_fasta.baseName}-SCAN
    """
    //
}
