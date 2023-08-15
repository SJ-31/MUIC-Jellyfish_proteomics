process PEPNET {
    publishDir "$outdir", mode: "copy"
    conda "/mnt/data/shannc/anaconda3/envs/PepNet"

    input:
    path(mgf)
    path(outdir)
    //

    output:
    path("${mgf.baseName}_pepnet.tsv")
    //

    script:
    """
    python $params.pepnet_exe \
        --input $mgf \
        --model $params.pepnetmodel \
        --output ${mgf.baseName}_pepnet.tsv
    """
    //
}

process EXTRACT_PEPNET {
    publishDir "$output", mode: "copy"

    input:
    path(pepnet_output)
    val(outdir)
    //

    output:
    path("pepnet*")
    //

    script:
    header="TITLE	DENOVO	Score	PPM Difference	Positional Score"
    """
    merge_tables.sh -r $header \
        -o all_pepnet.temp \
        -p tsv
    extract_denovo.py all_pepnet.temp pepnet
    """
    //
}
