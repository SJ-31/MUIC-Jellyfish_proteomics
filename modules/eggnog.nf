process EGGNOG {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/eggnog"

    input:
    path(unknown_fasta)
    val(outdir)
    //

    output:
    path("${unknown_fasta.baseName}*")
    //

    script:
    """
    export EGGNOG_DATA_DIR="$params.eggnog_data_dir"
    emapper.py -i $unknown_fasta \
        --output ${unknown_fasta.baseName} \
        --report_orthologs
    """
    //
}
