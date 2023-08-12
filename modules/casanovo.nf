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
        --config ${params.config}/casanovo.yaml
        --ouput=temp.tsv \
        --model=$params.casanovomodel
    grep -v ^M temp.tsv > ${mzMLs.baseName}_casanovo.tsv
    """
    //
}
