process CASANOVO {
    publishDir "$outdir", mode: "copy"
    conda "/mnt/data/shannc/anaconda3/casanovo"

    input:
    val(prefix)
    path(mzMLs)
    val(outdir)
    //
    output:
    path("${mzMLs.baseName}_casanovo.tsv")
    //

    script:
    """
    casanovo \
        --mode=denovo \
        --peak_path=$mzMLs \
        --ouput=${mzMLs.baseName}_casanovo.tsv \
        --model=$params.casanovomodel
    """
    //
}
