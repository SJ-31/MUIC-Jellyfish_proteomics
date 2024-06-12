process DEEPLOC {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/deeploc2"

    input:
    path(fasta)
    path(combined_results)
    path(intersected_searches)
    val(outdir)
    //

    output:
    path("deeploc_results.csv")
    path("${combined_results.baseName}.${combined_results.Extension}"), emit: tsv
    //

    script:
    def check = file("deeploc_results.csv")
    if (check.exists()) {
    """
    cp ${outdir}/deeploc_results.csv .
    """
    } else {
    """
    mkdir results
    deeploc2 -f $fasta \
            -o results
    cp results/* deeploc_results.csv

    if [[ -e deeploc_results.csv ]]; then
        Rscript $params.bin/R/merge_deeploc.r \
            -d deeploc_results.csv \
            -i $combined_results \
            -m $intersected_searches
    fi
    """
    }
    //
}
