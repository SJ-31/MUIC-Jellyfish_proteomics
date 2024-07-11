process PEPNET {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/PepNet"
    memory "10 GB"
    publishDir "$params.logs", mode: "copy", pattern: "*.log"

    input:
    path(mgf)
    val(outdir)
    //

    output:
    path(output), emit: peps
    path("*.log"), emit: log
    //

    script:
    output = "${mgf.baseName}_pepnet.tsv"
    check = file("${outdir}/${output}")
    if (check.exists()) {
        """
        cp "${outdir}/${output}" .
        cp "${outdir}/pepnet.log" .
        """
    } else {
    """
    python $params.pepnet_exe \
        --input $mgf \
        --model $params.pepnetmodel \
        --output ${mgf.baseName}_pepnet.tsv > pepnet.log
    """
    }
    //
}

process EXTRACT_PEPNET {
    publishDir "$outdir", mode: "copy"

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
    merge_tables.sh -r "$header" \
        -o all_pepnet.temp \
        -p tsv
    extract_denovo.py all_pepnet.temp pepnet
    """
    //
}
